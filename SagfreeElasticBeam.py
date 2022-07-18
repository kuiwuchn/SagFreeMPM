import time
import taichi as ti
import numpy as np
import sys
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import lsqr
from praxis import *

ti.init(arch=ti.x64) 

n_grid = 32
beam_x, beam_y = 6, 13
beam_width, beam_height = 20, 6

n_particles = beam_width * beam_height * 4

dx, inv_dx = 1 / float(n_grid), float(n_grid)
dt     = 1e-4 
gravity = 20.0
p_vol, p_rho = (dx * 0.5)**2, 1
p_mass = p_vol * p_rho
E, nu  = 0.1e4, 0.2                                                 # Young's modulus and Poisson's ratio
mu, la = E / (2 * (1 + nu)), E * nu / ((1+nu) * (1 - 2 * nu))       # Lame parameters
x      = ti.Vector.field(2,    dtype=float, shape=n_particles)      # position
v      = ti.Vector.field(2,    dtype=float, shape=n_particles)      # velocity
C      = ti.Matrix.field(2, 2, dtype=float, shape=n_particles)      # affine velocity field
F      = ti.Matrix.field(2, 2, dtype=float, shape=n_particles)      # deformation gradient
color  = ti.field(             dtype=int,   shape=n_particles)      # color id
Jp     = ti.field(             dtype=float, shape=n_particles)      # plastic deformation
grid_v = ti.Vector.field(2,    dtype=float, shape=(n_grid, n_grid)) # grid node momentum/velocity
grid_m = ti.field(             dtype=float, shape=(n_grid, n_grid)) # grid node mass
grid_i = ti.field(             dtype=int,   shape=(n_grid, n_grid)) # grid node index

####################################################################
# functions for the forward MPM simulation 

@ti.kernel
def cleanGrid():
    for i, j in grid_m:
        grid_v[i, j] = [0, 0]
        grid_m[i, j] = 0

@ti.kernel
def P2G():
    for p in x:
        base = (x[p] * inv_dx - 0.5).cast(int)
        fx = x[p] * inv_dx - base.cast(float)
        w = [0.5 * (1.5 - fx) ** 2, 0.75 - (fx - 1) ** 2, 0.5 * (fx - 0.5) ** 2]
        F[p] = (ti.Matrix.identity(ti.f32, 2) + dt * C[p]) @ F[p] 
        R, S = ti.polar_decompose(F[p])
        Jp[p] = Jp[p] * (1 + dt * C[p].trace()) 
        cauchy = (2 * mu * (F[p] - R) + la * (R.transpose() @ F[p] -  ti.Matrix.identity(ti.f32, 2)).trace() * R) @ F[p].transpose() 
        stress = (-dt * p_vol * 4 * inv_dx * inv_dx) * cauchy 
        affine = stress + p_mass * C[p]
        for i, j in ti.static(ti.ndrange(3, 3)): 
            offset = ti.Vector([i, j])
            dpos = (offset.cast(float) - fx) * dx
            weight = w[i][0] * w[j][1]
            grid_v[base + offset] += weight * (p_mass * v[p] + affine @ dpos)
            grid_m[base + offset] += weight * p_mass


@ti.kernel
def updateGrid(g : ti.f32):
    for i, j in grid_m:
        if grid_m[i, j] > 0:
            grid_v[i, j] = (1 / grid_m[i, j]) * grid_v[i, j] # Momentum to velocity
            grid_v[i, j][1] -= dt * g

        # boundary 
        if i < 3 and grid_v[i, j][0] < 0:          grid_v[i, j][0] = 0 
        if i > n_grid - 3 and grid_v[i, j][0] > 0: grid_v[i, j][0] = 0
        if j < 3 and grid_v[i, j][1] < 0:          grid_v[i, j][1] = 0
        if j > n_grid - 3 and grid_v[i, j][1] > 0: grid_v[i, j][1] = 0

        # fixed left side of the beam
        if i < beam_x + 3:  grid_v[i, j] = [0, 0]

@ti.kernel
def G2P():
    for p in x:
        base = (x[p] * inv_dx - 0.5).cast(int)
        fx = x[p] * inv_dx - base.cast(float)
        w = [0.5 * (1.5 - fx) ** 2, 0.75 - (fx - 1.0) ** 2, 0.5 * (fx - 0.5) ** 2]
        new_v = ti.Vector.zero(ti.f32, 2)
        new_C = ti.Matrix.zero(ti.f32, 2, 2)
        for i, j in ti.static(ti.ndrange(3, 3)): 
            dpos = ti.Vector([i, j]).cast(float) - fx
            g_v = grid_v[base + ti.Vector([i, j])]
            weight = w[i][0] * w[j][1]
            new_v += weight * g_v
            new_C += 4 * inv_dx * weight * g_v.outer_product(dpos)
        v[p], C[p] = new_v, new_C
        x[p] += dt * v[p]

####################################################################
# initialization function

@ti.kernel
def initialize():
    for i in range(beam_width * 2):
        for j in range(beam_height * 2):
            p_count = j + i * beam_height * 2
            x[p_count] = ti.Matrix([(float(beam_x) + i * 0.5 + 0.25) * dx, (float(beam_y) + j * 0.5 + 0.25) * dx])
            color[p_count] = 1
            v[p_count] = ti.Matrix([0, 0])
            F[p_count] = ti.Matrix([[1, 0], [0, 1]])
            Jp[p_count] = 1

####################################################################
# functions for sagfree initialization 

@ti.kernel
def P2GStatic():
    for p in x:
        base = (x[p] * inv_dx - 0.5).cast(int)
        fx = x[p] * inv_dx - base.cast(float)
        w = [0.5 * (1.5 - fx) ** 2, 0.75 - (fx - 1) ** 2, 0.5 * (fx - 0.5) ** 2]
        R, S = ti.polar_decompose(F[p])
        cauchy = (2 * mu * (F[p] - R) + la * (R.transpose() @ F[p] -  ti.Matrix.identity(ti.f32, 2)).trace() * R) @ F[p].transpose() 
        affine = (-p_vol * 4 * inv_dx * inv_dx) * cauchy
        for i, j in ti.static(ti.ndrange(3, 3)): 
            offset = ti.Vector([i, j])
            dpos = (offset.cast(float) - fx) * dx
            weight = w[i][0] * w[j][1]
            grid_v[base + offset] += weight * (affine @ dpos)
            grid_m[base + offset] += weight * p_mass

def computeRowAndColNum():
    for i in range(0, n_grid):
        for j in range(0, n_grid):
            grid_i[i, j] = -1

    rn = 0
    for i in range(0, n_grid):
        for j in range(0, n_grid):
            if grid_m[i, j] > 0 and i >= beam_x + 3:
                grid_i[i, j] = rn
                rn = rn + 1

    cn = 0
    for i in range(0, n_particles):
        if x[i][0] > (beam_x + 1) * dx:   # extra columns for particles
            cn = cn + 1   

    return rn, cn

def getSpMCoeff() :
    print("get sparse matrix information 2D")
    start_time = time.time()
    arr_len = 0
    for p in range(0, n_particles): 
        if x[p][0] > (beam_x + 1) * dx:
            basex = int(x[p][0] * inv_dx - 0.5)
            basey = int(x[p][1] * inv_dx - 0.5)
            fx = np.array([float(x[p][0] * inv_dx - basex), float(x[p][1] * inv_dx - basey)])
            w = [0.5 * (1.5 - fx) ** 2, 0.75 - (fx - 1) ** 2, 0.5 * (fx - 0.5) ** 2]
            for i in range(0, 3): 
                for j in range(0, 3): 
                    offset = np.array([i, j]).astype(np.float)
                    dpos = (offset - fx) * dx
                    weight = w[i][0] * w[j][1]
                    if (grid_m[basex + i, basey + j] > 0 and basex + i >= beam_x + 3) :
                        arr_len = arr_len + 1

    row = np.zeros(arr_len*4)
    col = np.zeros(arr_len*4)
    dat = np.zeros(arr_len*4)
    arr_len = 0  
    col_num = 0
    for p in range(0, n_particles): 
        if x[p][0] > (beam_x + 1) * dx:
            basex = int(x[p][0] * inv_dx - 0.5)
            basey = int(x[p][1] * inv_dx - 0.5)
            fx = np.array([float(x[p][0] * inv_dx - basex), float(x[p][1] * inv_dx - basey)])
            w = [0.5 * (1.5 - fx) ** 2, 0.75 - (fx - 1) ** 2, 0.5 * (fx - 0.5) ** 2]
            for i in range(0, 3): 
                for j in range(0, 3): 
                    offset = np.array([i, j]).astype(np.float)
                    dpos = (offset - fx) * dx
                    weight = w[i][0] * w[j][1]
                    if (grid_m[basex + i, basey + j] > 0 and basex + i >= beam_x + 3) :
                        row_id = grid_i[basex + i, basey + j]
                        col_id = col_num
                        row[arr_len*4+0]=row_id * 2 + 0
                        col[arr_len*4+0]=col_id * 3 + 0
                        dat[arr_len*4+0]=weight * dpos[0]
                   
                        row[arr_len*4+1]=row_id * 2 + 0
                        col[arr_len*4+1]=col_id * 3 + 2
                        dat[arr_len*4+1]=weight * dpos[1]
                     
                        row[arr_len*4+2]=row_id * 2 + 1
                        col[arr_len*4+2]=col_id * 3 + 2
                        dat[arr_len*4+2]=weight * dpos[0]
                         
                        row[arr_len*4+3]=row_id * 2 + 1
                        col[arr_len*4+3]=col_id * 3 + 1
                        dat[arr_len*4+3]=weight * dpos[1]
                        arr_len = arr_len + 1
            col_num = col_num + 1
    print("-- done")
    print("construction time: %s seconds" % (time.time() - start_time))
    return row, col, dat

def getRhsVec(row_num) :
    rhs = np.zeros(row_num * 2)
    row_num = 0
    for i in range(0, n_grid):
       for j in range(0, n_grid):
           if grid_m[i, j] > 0 and i >= beam_x + 3:
               rhs[row_num * 2 + 0] = 0
               rhs[row_num * 2 + 1] = gravity * grid_m[i, j]
               row_num = row_num + 1
    return rhs

def directSolverLsqr(A_sp, b) :
    print("direct solver - lsqr")
    start_time = time.time()
    c, istop, itn, normr = lsqr(A_sp, b)[:4]
    print("solver time: %s seconds" % (time.time() - start_time))
    return c

def evalAffine(f00, f11, f10):
    localF = np.array([[f00, f10], [f10, f11]])
    U, S, Vt = np.linalg.svd(localF)
    R = U @ Vt
    S = Vt.T @ np.diag(S) @ Vt
    cauchy = (2 * mu * (localF - R) + la * np.trace(R.T @ localF - np.eye(2)) * R) @ localF.T
    affine = (-p_vol * 4 * inv_dx * inv_dx) * cauchy
    return affine

def findSol(A0) :
    n = 3
    t0 = 1e-20
    h0 = 1
    prin = 0
    
    def f(r, n):
       A = evalAffine(r[0], r[1], r[2])                
       diff = A - A0
       val = diff[0, 0] * diff[0, 0]
       val += diff[0, 1] * diff[0, 1]
       val += diff[1, 0] * diff[1, 0]
       val += diff[1, 1] * diff[1, 1]
       return val

    r = np.array([1.0, 1.0, 0.0])
    pr, r = praxis ( t0, h0, n, prin, r, f )

    return r

####################################################################
# functions for verification 

@ti.kernel
def P2GStaticTest():
    for p in x:
        base = (x[p] * inv_dx - 0.5).cast(int)
        fx = x[p] * inv_dx - base.cast(float)
        w = [0.5 * (1.5 - fx) ** 2, 0.75 - (fx - 1) ** 2, 0.5 * (fx - 0.5) ** 2]
        affine = F[p]
        for i, j in ti.static(ti.ndrange(3, 3)): 
            offset = ti.Vector([i, j])
            dpos = (offset.cast(float) - fx) * dx
            weight = w[i][0] * w[j][1]
            grid_v[base + offset] += weight * (affine @ dpos)
            grid_m[base + offset] += weight * p_mass

def verifyGlobalStep():
    cleanGrid()
    col_num = 0
    for p in range(0, n_particles): 
        if x[p][0] >= (beam_x + 1) * dx:
            F[p][0, 0] = c[col_num * 3 + 0]
            F[p][1, 1] = c[col_num * 3 + 1]
            F[p][0, 1] = c[col_num * 3 + 2]
            F[p][1, 0] = c[col_num * 3 + 2]
            col_num = col_num + 1
    
    P2GStaticTest()

    for i in range(0, n_grid):
        for j in range(0, n_grid):
            if grid_m[i, j] > 0 and i >= beam_x + 3:
                if not np.isclose(grid_v[i, j][0], 0, atol = 1e-6) :
                    print("Error in global step: ", i, j, grid_v[i, j][0] - 0)
                if not np.isclose(grid_v[i, j][1], gravity * grid_m[i, j], atol = 1e-6) :
                    print("Error in global step: ", i, j, grid_v[i, j][1] - gravity * grid_m[i, j])

def verifyLocalStep(res_array, col_num):
    for i in range(0, col_num):  
        affine = evalAffine(res_array[i * 3 + 0], res_array[i * 3 + 1], res_array[i * 3 + 2])
        res_vec = np.array([affine[0, 0], affine[1, 1], affine[1, 0]])
        tar_vec = np.array([c[i * 3], c[i * 3 + 1], c[i * 3 + 2]])
        if not np.allclose(res_vec, tar_vec, atol = 1e-6):
            print("Error in local step: ", c[i * 3], c[i * 3 + 1], c[i * 3 + 2], res_vec, tar_vec, res_vec - tar_vec)
            
####################################################################
# start from here
####################################################################

initialize()

# get lhs matrix
P2GStatic()
row_num, col_num = computeRowAndColNum()
row, col, dat = getSpMCoeff()
A_sp  = csc_matrix((dat, (row, col)), shape=(row_num * 2, col_num * 3))
At_sp = csc_matrix((dat, (col, row)), shape=(col_num * 3, row_num * 2))

# get rhs vector
b = getRhsVec(row_num)

# solve the global stage
c = directSolverLsqr(A_sp, b)      

# verify results from the global stage
verifyGlobalStep()

# solve the local stage
res_array = np.zeros(col_num * 3)
for i in range(0, col_num) :
    res = findSol(np.array([[c[i * 3], c[i * 3 + 2]], [c[i * 3 + 2], c[i * 3 + 1]]]))
    res_array[i * 3 + 0] = res[0]
    res_array[i * 3 + 1] = res[1]
    res_array[i * 3 + 2] = res[2]
  
# verify results from the local stage  
verifyLocalStep(res_array, col_num)

# copy the results into F
col_num = 0
for p in range(0, n_particles): 
    if x[p][0] >= (beam_x + 1) * dx:
        F[p][0, 0] = res_array[col_num * 3 + 0]
        F[p][1, 1] = res_array[col_num * 3 + 1]
        F[p][0, 1] = res_array[col_num * 3 + 2]
        F[p][1, 0] = res_array[col_num * 3 + 2]
        col_num = col_num + 1

gui = ti.GUI("Sagfree elastic beam", res=512, background_color=0x222222)

frame = 0
while not gui.get_event(ti.GUI.ESCAPE, ti.GUI.EXIT):

    print("frame ", frame)
    if frame > 200 : 
        gravity = -20.0 * np.sin(frame)
    frame = frame + 1

    for s in range(int(2e-3 // dt)):
        cleanGrid()
        P2G()
        updateGrid(gravity)
        G2P()
    gui.circles(x.to_numpy(), radius=1.5, color=0xED553B)

    gui.show()

    cleanGrid()
