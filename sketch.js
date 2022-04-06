//source https://mikeash.com/pyblog/fluid-simulation-for-dummies.html

const N = 100;
const iter = 4;
const SCALE = 5; 

let fluid;

function setup() {
  createCanvas(N*SCALE,N*SCALE);
  //frameRate(30);
  fluid = new Fluid(1e-1,0,0);
}

function mouseDragged(){
  fluid.addDensity(floor(mouseX/SCALE)%N, floor(mouseY/SCALE)%N, 1e-8);
  //fluid.addVelocity(floor(mouseX/SCALE)%N, floor(mouseY/SCALE)%N, 1e-8, 1e-8);
  //console.log("clock");
}

function draw() {
  background(0);
  fluid.step();
  fluid.renderD();
}
