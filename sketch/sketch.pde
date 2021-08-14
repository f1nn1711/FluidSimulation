int size = 125;
int cellSize = 5;
int itterations = 10;

Fluid fluid;

int nSources = 5;

boolean sourceOn = true;

int view = 0;
int maxView = 2;

ArrayList<Source> sources = new ArrayList<Source>();

void settings() {
  size(size*cellSize,size*cellSize);
}

void setup() {
  fluid = new Fluid();
  
  for (int n = 1; n < nSources; n++) {
    int y = (floor(size/(nSources)))*n;
    sources.add(
      new Source(5, y, 6, 0, 500)
    );
  }
}

void draw() {
  background(0);
  
  fluid.step();
  fluid.dampen();
  fluid.render();
  
  if (sourceOn == true){
    for (Source s : sources) {
      fluid.addConc(s.x, s.y, s.qty);
      fluid.addVel(s.x, s.y, s.xComp, s.yComp);
    }
  }
  

}

void mouseDragged() {
  //fluid.addConc(mouseX/cellSize, mouseY/cellSize, 100);
  fluid.addVel(mouseX/cellSize, mouseY/cellSize, 2*(mouseX-pmouseX), 2*(mouseY-pmouseY));
  
  //fluid.addConc(20, height/2, 100);
  //fluid.addVel(20, height/2, 30, 1);
}

void keyPressed() {
  switch (keyCode) {
    case 9://tab
      if (sourceOn == true) {
        sourceOn = false;
      } else {
        sourceOn = true;
      }
      break;
    case 86://v
      if (view != maxView) {
        view++;
      } else {
        view = 0;
      }
  }

}
