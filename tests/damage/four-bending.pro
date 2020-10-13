
mpart.overlap = 1; // for parallel simulations

log =
{
  // Print informational messages to the terminal.
  pattern = "*.info";
  file    = "-$(CASE_NAME).log";
};


input =
{
  // Specify the name of the input data file.
  file = "four-bending.data";
};

control =
{
  pause = 0.;
  runWhile = "i < 0";
};


model =  "Matrix"
{
  matrix.type = "Sparse";
  model       = "Multi"
  {
    models = [ "all", "load", "lodi", "monitor" ];

    all = "NonLocalDamage"
    {
      elements  = "bulk";
      thickness = 50.;

      // Material parameters.

      materials = ["mat1"];
      mat1 = 
      {
	       type       = "Damage";
           dim        = 2;
	       state      = "PLANE_STRESS";   
           young      = 40e3;
           poisson    = 0.2;
           softening  = "exponential1";
	       equistrain = "vonMises";
	       kappaI     = 0.000075;
           alpha      = 0.92;
           beta       = 300.;
           eta        = 10.;
      
           elements  = "bulk";
      };
      
      // The shape object to be used. This depends on the finite
      // element mesh.

      shape =
      {
	    type        = "Triangle3";
        intScheme   = "Gauss1";
      };

      radius = 2.8284;
      weightFunc = " let radius = 2.8284;"
            " return (1.0 - (r/radius)^2)^2"; averaging = "symmetric";
    };
 
    load =  "LoadScale"
    {
        scaleFunc = "1e-3 * (i)";
  
        model =
        {
          type     = "Constraints";
          conTable = "load";
        };
    };

    // responsible for load-displacement curves

    lodi =
    {
      type   = "Lodi";
      group  = "disp";
    };    

    monitor = {
      type = "Monitor";
      node = 17; 
    };
  };
};
 
extraModules =
{
  modules = ["solver", "graph", "lodi", "vtk", "view"];

  solver =  "Nonlin"
  {
     precision = 1.0e-5;
     maxIter   = 20;
     solver  = "SparseLU"
     {
         useThreads = true;
     };
  };

  graph =
  {
      type = "Graph";
      dataSets = "loadDisp";
      loadDisp =
      {
        key = "Load-displacement curve";
        xData = "-model.model.monitor.state.dy"; 
        yData = "-0.05*model.model.lodi.load[1]";
      };
   };
    
   // if you want to write the load-displacement data to file
    lodi = "Sample"
    {
      file       = "four-bending-lodi1.dat";
      header     = "  uy | fy ";
      dataSets   = ["-model.model.monitor.state.dy","-0.05* model.model.lodi.load[1]"];
    };
    
    vtk = "Paraview"
    {
       fileName   = "four-bending";
       elements = "bulk";
       interval = 50;
       dofs     = ["dx","dy"];
       data     = ["damage"];
    };

view = "FemView"
{
  window =
  {
     height = 300;
     width  = 900;
  };
  
  //snapFile = "$(CASE_NAME)%2i.png";
  configFile  = "$(CASE_NAME).view";
  //snapWhen = "i%20";

  // Define some data sets.

  dataSets = [ "disp" , "damage" ];

  disp =
  {
    type   = "Vector";
    vector = "state";
   };

   damage =
   {
      type  = "Table";
      table = "nodes/damage";
   };

   mesh =
   {
       deformation = "80 * disp";
       plugins = "colors";
       colors = {
          type = "MeshColorView";
          data = "damage";
          palette   = "custom";
          autoScale = false;
      };
    };
};
};
