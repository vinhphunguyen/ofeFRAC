
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
  file = "$(CASE_NAME).data";
};

control =
{
  pause = 0.;
  runWhile = "i < 0";
};


model =  "Matrix"
{
  // The outer model is responsible for creating the global
  // stiffness matrix.

  matrix.type = "FEM";

  model       = "Multi"
  {
    // The inner models is responsible for assembling the global
    // stiffness matrix. This model is a so-called MultiModel.

    models = [ "all", "load", "lodi"];

    all = "GradientDamage"
    {
      elements  = "bulk";
      thickness = 50.;

      // Material parameters.

      material = 
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
           beta       = 500.;
           eta        = 10.;
           lengthscale= 1.0;
      };
      
      // The shape object to be used. This depends on the finite
      // element mesh.

      shape =
      {
	    type        = "Triangle";
        intScheme   = "Gauss1";
      };

      cMatrix =
      {
        type = "StressBased";
        //type = "Constant";
        c0   = 1.;
        ft   = 3.;
      };
    };

     
    load =
    {
      type = "EnergyArclen";

      arcFunc = "ERC";

      swtIter = 3;
      optIter = 5;

      // The initial load increment:

      loadIncr = 1.0;

      // The maximum and minimum absolute load increments:

      minIncr = 0.0005;
      maxIncr = 0.5;

      // The initial load scale:

      loadScale = 0.0;

      model =
      {
	      type      = "PointLoad";
	      loadTable = "load";
      };
    };

    // responsible for load-displacement curves

    lodi =
    {
      type    = "LoadCMOD";
      fNodes  = "disp";
      cmoNodes  = [17,2];
      file    = "$(CASE_NAME).dat";
    };    

  };
};
 
extraModules =
{
  modules = ["solver", "graph", "lodi", "vtk"];
  
    solver = "Nonli"  // can be Nonlin, Linsolve...
    {
      precision = 1.0e-6;
    
      maxIter = 20;

      solver =
      {
        type = "GMRES"; 
        //type="SparseLU";
        precon.type="ILUd";
        useThreads=true;
        precision = 1.0e-6;
      };
    };


  graph =
  {
      type = "Graph";
      dataSets = "loadDisp";
      loadDisp =
      {
        key = "Load-displacement curve";
        xData = "-model.model.lodi.disp[1]"; 
        yData = "-0.05*model.model.lodi.load[1]";
      };
   };
    
   vtk = "Paraview"
   {
      fileName = "kalthoff";
      elements = "bulk";
      data = ["stress","strain"]; 
      interval = 50;
  };
}; // end of extraModules

view =
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
