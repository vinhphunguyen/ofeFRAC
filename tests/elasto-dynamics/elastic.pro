log =
{
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};

control = 
{
  pause = 0.;
  runWhile = "i < 0";
};

input =
{
  file = "cylinder.data";
};

model = "Matrix"
{
  matrix.type = "FEM";

  model       =  "Multi"
  {
    models = [ "bulk",  "load" ];

    bulk = "Elasticity"
    {
      elements = "bulk";

      shape =
      {
        type      = "Tet4";
        intScheme = "Gauss1";
      };
      
      material = ["ddd"];
      
      ddd =
      {
        type   = "Hooke";
        dim    = 3;
        young      = 210e3;
        poisson    = 0.2;
	      rho   = 7850e-12;

	      elements   = "bulk";
      };
    };

      load = 
        {
          type = "LoadScale";
    
          scaleFunc = "return
               exp(-time/1e-4)";

	       //scaleFunc = "1e-2 * (i)";
    
          model = "LineLoad"
          {
            elements = "inner";
            load     = 1.;
            shape=
            {
              type = "BTriangle3";
            };
          };
        };   

  };
};


extraModules =
{
  modules = ["solver","vtk"];

// Implicit dynamics
  
solver= "Newmark" 
{
  deltaTime = 1e-5; //initial time step

  solver = "Nonlin" 
  {
    precision = 1.0e-4;
    maxIter   = 100;
  
    solver =
    {
      type="AGMRES";
      precision = 1.0e-7;

      // This is for some useful solver information.
      //noiseLevel=1;

      // The larger this value, the more memory is used
      // but, in general, the convergence rate is improved.
      restartIter=150;

      precon = "ILUd"
      {
        // Determines the memory size of the preconditioner with
        // respect to the size of the global stiffness matrix.
        // The larger this value, the more effective the
        // preconditioner (setting it to a very large value
        // turns the preconditioner into a direct solver).
        maxFill = 5.0;

        // A larger value increases the cost of setting up the
        // preconditioner, but also makes the preconditioner
        // more effective.
        quality = 2.5;
      };
    };
  };
};


  vtk =  "Paraview"
  {
      fileName   = "$(CASE_NAME)";
      elements = "bulk";
      interval = 1;
      data     = ["stress"];
      dofs     = ["dx","dy"];
  };
};


