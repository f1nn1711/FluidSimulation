class Fluid {
  float timeStep = 0.01;
  float diffusion = 0.000001;
  float viscosity = 0.0000001;
  float dampening = 0.000000000001;

  float[] prevDensity;
  float[] density;

  float[] prevVx;
  float[] prevVy;

  float[] Vx;
  float[] Vy;

  Fluid() {
    this.prevDensity = new float [size*size];
    this.density = new float [size*size];

    this.prevVx = new float [size*size];
    this.prevVy = new float [size*size];

    this.Vx = new float [size*size];
    this.Vy = new float [size*size];
  }

  void addConc(int x, int y, float qty) {
    int ref = getRef(x, y);
    this.density[ref] += qty;
  }

  void addVel(int x, int y, float xComp, float yComp) {
    int ref = getRef(x, y);
    this.Vx[ref] += xComp;
    this.Vy[ref] += yComp;
  }
  
  void step() {
    float[] Vx      = this.Vx;
    float[] Vy      = this.Vy;
    
    float[] Vx0     = this.prevVx;
    float[] Vy0     = this.prevVy;

    float[] s       = this.prevDensity;
    float[] density = this.density;
    
    diffuse(1, Vx0, Vx, this.viscosity, this.timeStep);
    diffuse(2, Vy0, Vy, this.viscosity, this.timeStep);
    
    project(Vx0, Vy0, Vx, Vy);
    
    advect(1, Vx, Vx0, Vx0, Vy0, this.timeStep);
    advect(2, Vy, Vy0, Vx0, Vy0, this.timeStep);
    
    project(Vx, Vy, Vx0, Vy0);
    
    diffuse(0, s, density, this.diffusion, this.timeStep);
    advect(0, density, s, Vx, Vy, this.timeStep);  
  }
  
  void render() {
    switch (view) {
      case 0:
        for (int i = 0; i < size; i++) {
          for (int j = 0; j < size; j++) {
            float x = i * cellSize;
            float y = j * cellSize;
            float d = this.density[getRef(i,j)];
            
            fill(d);
            noStroke();
            square(x, y, cellSize);
          }
        }
        break;
      case 1:
        float[] mapVx = new float[this.Vx.length];
        float[] mapVy = new float[this.Vy.length];
        
        arrayCopy(this.Vx, mapVx);
        arrayCopy(this.Vy, mapVy);
        
        float minVx = min(mapVx);
        float minVy = min(mapVy);
        
        float maxVx = max(mapVx);
        float maxVy = max(mapVy);
        
        int everyOtherX = 2;
        int everyOtherY = 2;
        
        for (int i = 0; i < mapVx.length; i++) {
          mapVx[i] = map(mapVx[i], minVx, maxVx, 0.1, cellSize*everyOtherX);
          mapVy[i] = map(mapVy[i], minVy, maxVy, 0.1, cellSize*everyOtherY);
        }
      
        for (int i = 0; i < size; i += everyOtherX) {
          for (int j = 0; j < size; j += everyOtherY) {
            float x = (i * cellSize)+(cellSize/2);
            float y = (j * cellSize)+(cellSize/2);

            float xComp = mapVx[getRef(i,j)];
            float yComp = mapVy[getRef(i,j)];
            stroke(255, 0, 0);
            strokeWeight(1);
            line(x, y, (x+xComp), (y+yComp));
          }
        }
        break;
      case 2:
        float[] dCopy = new float[this.density.length];
        arrayCopy(this.density, dCopy);
        
        float maxD = max(dCopy);
        float minD = min(dCopy);
        
        for (int dn = 0; dn < dCopy.length; dn++) {
          dCopy[dn] = map(dCopy[dn], minD, maxD, 0, 1);

        }
      
        for (int i = 0; i < size; i++) {
          for (int j = 0; j < size; j++) {
            float x = i * cellSize;
            float y = j * cellSize;
            
            float d = dCopy[getRef(i,j)];
            
            //int[] rgbValues = greyToRGB(d);
            int[] rgbValues = HeatColor(d, minD, maxD);
            fill(rgbValues[0], rgbValues[1], rgbValues[2]);
            noStroke();
            square(x, y, cellSize);
          }
        }
        break;
    }
  }
  
  void dampen() {
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < size; j++) {
        this.density[getRef(i,j)] *= 1-this.dampening;
      }
    }
  }
}

public void diffuse(int b, float[] x, float[] prevX, float diffusion, float timestep) {
  float a = timestep * diffusion * (size - 2) * (size - 2);
  linearSolution(b, x, prevX, a, 1 + 6 * a);
};

public void linearSolution(int b, float[] x, float[] x0, float a, float c)
{
  float cRecip = 1.0 / c;
  for (int k = 0; k < itterations; k++) {
    for (int j = 1; j < size - 1; j++) {
      for (int i = 1; i < size - 1; i++) {
        x[getRef(i, j)] =
          (x0[getRef(i, j)]
          + a*(    x[getRef(i+1, j)]
          +x[getRef(i-1, j)]
          +x[getRef(i, j+1)]
          +x[getRef(i, j-1)]
          )) * cRecip;
      }
    }
    set_bnd(b, x);
  }
}

public void project(float[] Vx, float[] Vy, float[] p, float[] div) {
  for (int j = 1; j < size - 1; j++) {
    for (int i = 1; i <size - 1; i++) {
      div[getRef(i, j)] = -0.5f*(
        Vx[getRef(i+1, j)]
        -Vx[getRef(i-1, j)]
        +Vy[getRef(i, j+1)]
        -Vy[getRef(i, j-1)]
        )/size;
      p[getRef(i, j)] = 0;
    }
  }

  set_bnd(0, div); 
  set_bnd(0, p);
  linearSolution(0, p, div, 1, 6);

  for (int j = 1; j < size - 1; j++) {
    for (int i = 1; i < size - 1; i++) {
      Vx[getRef(i, j)] -= 0.5f * (p[getRef(i+1, j)] - p[getRef(i-1, j)]) * size;
      Vy[getRef(i, j)] -= 0.5f * (p[getRef(i, j+1)] - p[getRef(i, j-1)]) * size;
    }
  }
  set_bnd(1, Vx);
  set_bnd(2, Vy);
}

public void advect(int b, float[] d, float[] prevD, float[] Vx, float[] Vy, float timestep)
{
  float i0, i1, j0, j1;

  float dtx = timestep * (size - 2);
  float dty = timestep * (size - 2);

  float s0, s1, t0, t1;
  float tmp1, tmp2, x, y;

  float ifloat, jfloat;
  int i, j;

  for (j = 1, jfloat = 1; j < size - 1; j++, jfloat++) { 
    for (i = 1, ifloat = 1; i < size - 1; i++, ifloat++) {
      tmp1 = dtx * Vx[getRef(i, j)];
      tmp2 = dty * Vy[getRef(i, j)];

      x = ifloat - tmp1; 
      y = jfloat - tmp2;


      if (x < 0.5f) x = 0.5f; 
      if (x > size + 0.5f) x = size + 0.5f; 
      i0 = floor(x); 
      i1 = i0 + 1.0f;
      if (y < 0.5f) y = 0.5f; 
      if (y > size + 0.5f) y = size + 0.5f; 
      j0 = floor(y);
      j1 = j0 + 1.0f; 

      s1 = x - i0; 
      s0 = 1.0f - s1; 
      t1 = y - j0; 
      t0 = 1.0f - t1;


      int i0i = int(i0);
      int i1i = int(i1);
      int j0i = int(j0);
      int j1i = int(j1);


      d[getRef(i, j)] = 
        s0 * ( t0 * prevD[getRef(i0i, j0i)] + t1 * prevD[getRef(i0i, j1i)]) +
        s1 * ( t0 * prevD[getRef(i1i, j0i)] + t1 * prevD[getRef(i1i, j1i)]);
    }
  }

  set_bnd(b, d);
}

public void set_bnd(int b, float[] x)
{
  for (int i = 1; i < size - 1; i++) {
    x[getRef(i, 0)] = b == 2 ? -x[getRef(i, 1)] : x[getRef(i, 1)];
    x[getRef(i, size-1)] = b == 2 ? -x[getRef(i, size-2)] : x[getRef(i, size-2)];
  }

  for (int j = 1; j < size - 1; j++) {
    x[getRef(0, j)] = b == 1 ? -x[getRef(1, j)] : x[getRef(1, j)];
    x[getRef(size-1, j)] = b == 1 ? -x[getRef(size-2, j)] : x[getRef(size-2, j)];
  }

  x[getRef(0, 0)] = 0.5f * (x[getRef(1, 0)] + x[getRef(0, 1)]);
  x[getRef(0, size-1)] = 0.5f * (x[getRef(1, size-1)] + x[getRef(0, size-2)]);
  x[getRef(size-1, 0)] = 0.5f * (x[getRef(size-2, 0)] + x[getRef(size-1, 1)]);
  x[getRef(size-1, size-1)] = 0.5f * (x[getRef(size-2, size-1)] + x[getRef(size-1, size-2)]);
}
