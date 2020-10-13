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
  file = "$(CASE_NAME).data";
};

model = "Matrix"
{
  matrix.type = "FEM";

  model       =  "Multi"
  {
    models = [ "bulk",  "force" ];

    bulk = "Elasticity"
    {
      elements = "bulk";

      shape =
      {
        type      = "Quad4";
        intScheme = "Gauss2*Gauss2";
      };
      
      material = ["ddd"];
      
      ddd =
      {
        type   = "Hooke";
        dim    = 2;
        state  = "PLANE_STRESS";

        young      = 1e3;
        poisson    = 0.2;

	    elements   = "bulk";
      };
    };
    

    force =  "LoadScale"
    {
       // time in seconds
       //scaleFunc = "exp(-time/1e-4)";

       model =
       {
         type     = "LineLoad";
         elements = "rightEdge";
         load     = 10.;
         shape=
         {
           type = "BLine2";
         };
       };
    };
  };
};

extraModules =
{
  modules = ["solver", "paraview", "view"];
  
  solver = 
  {
      type = "Nonlin";
    
      precision = 1.0e-5;
    
      maxIter   = 20;
  
      solver =
      {
        type = "SparseLU"; 
      };
  };

  paraview =  "Paraview"
  {
      fileName   = "plate-hole";
      elements = "bulk";
      interval = 1;
      data     = ["ss"];
      dofs=["dx","dy"];
  };

view = "FemView"
{
 
  window =
  {
    height = 300;
    width  = 900;
  };

  snapFile = "$(CASE_NAME)%2i.png";
  configFile  = "$(CASE_NAME).view";


  // Define some data sets.

  dataSets = [ "disp" , "stress" ];

  disp =
  {
    type   = "Vector";
    vector = "state";
  };

  stress =
  {
    type  = "Table";
    table = "nodes/stress";
  };

  mesh =
  {
    // Use the data set "disp" to deform the mesh.

    deformation = "1 * disp";

    // Define extra plugins for showing more data sets.

    plugins = "colors";

    colors =
    {
      type = "MeshColorView";
      data = "disp";
      palette   = "custom";
      autoScale = true;
    };
  };

};
};
