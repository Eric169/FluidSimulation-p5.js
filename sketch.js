//source https://mikeash.com/pyblog/fluid-simulation-for-dummies.html

const N = 64;
const iter = 4;
const SCALE = 10; 

class Fluid{
  constructor(dt, diffusion, viscosity){
    this.dt = dt;
    this.diff = diffusion;
    this.visc = viscosity;

    this.s = [];
    this.s.length = N*N;
    this.s.fill(0);

    this.density = [];
    this.density.length = N*N;
    this.density.fill(0);

    for(let i = 0; i<N*N; i++) console.log(this.density[i]);

    this.Vx = [];
    this.Vx.length = N*N;
    this.Vx.fill(0);
    this.Vy = [];
    this.Vy.length = N*N;
    this.Vy.fill(0);

    this.Vx0 = [];
    this.Vx0.length = N*N;
    this.Vx0.fill(0);
    this.Vy0 = [];
    this.Vy0.length = N*N;
    this.Vy0.fill(0);
  }

  IX(x,y){
   return x + y * N;
  }
  step(){
    let visc     = this.visc;
    let diff     = this.diff;
    let dt       = this.dt;
    let Vx      = this.Vx;
    let Vy      = this.Vy;
    let Vx0     = this.Vx0;
    let Vy0     = this.Vy0;
    let s       = this.s;
    let density = this.density;
    
    this.diffuse(1, Vx0, Vx, visc, dt);
    this.diffuse(2, Vy0, Vy, visc, dt);
    
    this.project(Vx0, Vy0, Vx, Vy);
    
    this.advect(1, Vx, Vx0, Vx0, Vy0, dt);
    this.advect(2, Vy, Vy0, Vx0, Vy0, dt);
    
    this.project(Vx, Vy, Vx0, Vy0);
    
    this.diffuse(0, s, density, diff, dt);
    this.advect(0, density, s, Vx, Vy);
  }

  addDensity(x,y,amount){
    let idx = this.IX(x,y);
    this.density[idx] += amount;
  }

  addVelocity(x,y,amountX, amountY){
    let idx = this.IX(x,y);
    this.Vx[idx] += amountX;
    this.Vy[idx] += amountY;
  }

  diffuse(b, x, x0, diff, dt){ //knows how to diffuse every x knowing the diffusion time and iterations
    let a = dt * diff * (N-2) * (N-2);
    this.lin_solve(b,x,x0,a,1 + 6*a);
  }

  lin_solve(b, x, x0, a, c){ //solves e linear equation
    //for every x it checks all the neighbours
    let cRecip = 1.0 / c;
    for(let j = 1; j < N-1; j++){
      for(let i = 1; i < N-1; i++){
        x[this.IX(i,j)] = 
          (x0[this.IX(i,j)]
            + a*( x[this.IX(i+1,j)]
            +x[this.IX(i-1,j)]
            +x[this.IX(i,j+1)]
            +x[this.IX(i,j-1)]
            )) * cRecip;
      }
    }

    this.set_bnd(b,x);
  }

  set_bnd(b, x){
    for(let k = 1; k < N - 1; k++) {
        for(let i = 1; i < N - 1; i++) {
            x[this.IX(i, 0)] = b == 2 ? -x[this.IX(i, 1)] : x[this.IX(i, 1)];
            x[this.IX(i, N-1)] = b == 2 ? -x[this.IX(i, N-2)] : x[this.IX(i, N-2)];
        }
    }
    for(let k = 1; k < N - 1; k++) {
        for(let j = 1; j < N - 1; j++) {
            x[this.IX(0  , j)] = b == 1 ? -x[this.IX(1  , j)] : x[this.IX(1  , j)];
            x[this.IX(N-1, j)] = b == 1 ? -x[this.IX(N-2, j)] : x[this.IX(N-2, j)];
        }
    }
    
    x[this.IX(0, 0)] = 0.5 * (x[this.IX(1, 0)] + x[this.IX(0, 1)]);
    x[this.IX(0, N-1)] = 0.5 * (x[this.IX(1, N-1)] + x[this.IX(0, N-2)]);
    x[this.IX(N-1, 0)] = 0.5 * (x[this.IX(N-2, 0)] + x[this.IX(N-1, 1)]);
    x[this.IX(N-1, N-1)]   = 0.5 * (x[this.IX(N-2, N-1)] + x[this.IX(N-1, N-2)]);
  }

  project(velocX, velocY, p, div){
    for(let j = 1; j<N-1; j++){
      for(let i = 1; i<N-1; i++){
        div[this.IX(i,j)] = -.5*(
            velocX[this.IX(i+1,j)]
          -velocX[this.IX(i-1,j)]
          +velocY[this.IX(i,j+1)]
          -velocY[this.IX(i,j-1)]
        )/N;
        p[this.IX(i,j)] = 0;
      }
    }
    this.set_bnd(0,div);
    this.set_bnd(0,p);
    this.lin_solve(0,p,div,1,6);

    for(let j = 1; j<N-1; j++){
      for(let i = 1; i<N-1; i++){
        velocX[this.IX(i,j)] -= .5 * (p[this.IX(i+1,j)] - p[this.IX(i-1,j)]) * N;
        velocY[this.IX(i,j)] -= .5 * (p[this.IX(i,j+1)] - p[this.IX(i,j-1)]) * N;
      }
    }

    this.set_bnd(1,velocX);
    this.set_bnd(2,velocY);
  }

  advect(b, d, d0, velocX, velocY, dt)
  {
    let i0, i1, j0, j1;
    
    let dtx = dt * (N - 2);
    let dty = dt * (N - 2);
    
    let s0, s1, t0, t1;
    let tmp1, tmp2, x, y;
    
    let Nfloat = N;
    let ifloat, jfloat;
    let i, j;
    
    for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
        for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            tmp1 = dtx * velocX[this.IX(i, j)];
            tmp2 = dty * velocY[this.IX(i, j)];
            x    = ifloat - tmp1; 
            y    = jfloat - tmp2;
            
            if(x < 0.5) x = 0.5; 
            if(x > Nfloat + 0.5) x = Nfloat + 0.5; 
            i0 = floor(x); 
            i1 = i0 + 1.0;
            if(y < 0.5) y = 0.5; 
            if(y > Nfloat + 0.5) y = Nfloat + 0.5; 
            j0 = floor(y);
            j1 = j0 + 1.0; 
            
            s1 = x - i0; 
            s0 = 1.0 - s1; 
            t1 = y - j0; 
            t0 = 1.0 - t1;
            
            let i0i = i0;
            let i1i = i1;
            let j0i = j0;
            let j1i = j1;
            
            d[this.IX(i, j)] = 
                s0 * ( t0 * d0[this.IX(i0i, j0i)] ) +( t1 * d0[this.IX(i0i, j1i)]) +
                s1 * ( t0 * d0[this.IX(i1i, j0i)]) +( t1 * d0[this.IX(i1i, j1i)]);
        }
    }
    this.set_bnd(b, d);
  }

  renderD(){  
    for(let i = 0; i<N; i++){
      for(let j = 0; j<N; j++){
        let x = i * SCALE;
        let y = j * SCALE;
        let d = this.density[this.IX(i,j)];
        //console.log(d);
        fill(d);
        noStroke();
        rect(x,y,SCALE,SCALE);
      }
    }
  }
}

var fluid = new Fluid(0.1,0,0);

function setup() {
  createCanvas(N*SCALE,N*SCALE);
}

function mouseDragged(){
  fluid.addDensity(mouseX/SCALE, mouseY/SCALE, 1);
}

function draw() {
  background(0);
  fluid.step();
  fluid.renderD();
}
