class Source {
  int x, y, xComp, yComp;
  float qty;
  
  Source(int sx, int sy, int xComponent, int yComponent, float amount) {
    x = sx;
    y = sy;
    xComp = xComponent;
    yComp = yComponent;
    qty = amount;
  }
}

public int getRef(int x, int y) {
  x = constrain(x, 0, size-1);
  y = constrain(y, 0, size-1);
  return x + y * size;
}

public int[] greyToRGB(float gsValue) {
  //gsValue will be between 0 and 1
  int rValue = round(255*gsValue);
  int gValue = round(255*(1-gsValue));
  int bValue = 0;
  return new int[] {rValue, gValue, bValue};
}

public int[] HeatColor(float v,float vmin,float vmax){
   int r=0, g=0, b=0;
   float x = (v-vmin)/(vmax-vmin);
   r = round(255*constrain(-4*abs(x-0.75) + 1.5,0,1));
   g = round(255*constrain(-4*abs(x-0.50) + 1.5,0,1));
   b = round(255*constrain(-4*abs(x-0.25) + 1.5,0,1));
   return new int[] {r, g, b};
}

public int[][] calculateLinePoints(int[][] point){
  int[][] points;
  
  points[0] = new int[] {1,4};
  
  return points;
}
