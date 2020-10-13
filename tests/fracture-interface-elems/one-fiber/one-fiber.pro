mpart.overlap=1;

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

model =  "MP"
{
model =  "Matrix"
{
  matrix.type = "Sparse";
  model       = "Multi"
  {
    models = [ "bulk",  "interface","disp",  "lodi" ];

    bulk = "Elasticity"
    {
      elements = "bulk";

      shape =
      {
        type      = "Triangle3";
        intScheme = "Gauss1";
      };
      
      material = ["matrix","fiber"];
      
      matrix =
      {
        type   = "Hooke";
        dim    = 2;
        state  = "PLANE_STRAIN";

        young      = 4e3;
        poisson    = 0.4;

	      elements   = "matrix";
      };

      fiber =
      {
        type   = "Hooke";
        dim    = 2;
        state  = "PLANE_STRAIN";

        young      = 4e4;
        poisson    = 0.33;

     	elements   = "fiber";
      };
    };

    interface =
    {
       type     = "Interface";
       meshFile = "one-fiber-interface.mesh";
       elements = "interfaces";
       shape =
       {
          type      = "BLine2";
          intScheme = "Gauss2"   ; // NC2
       };

       //dNodeGroup   = "right";
       //dirNodeGroup = "left";
       //dirNodeGroup = "node12";

       material = ["mat0"];
       mat0={
          type  = "TuronMixed";
	      dim   = 2;
	      dummy = 1e6;
	      gI    = 0.05;
	      gII   = 0.05;
	      f2t   = 10.;
	      f6    = 10.;
	      eta   = 2.;

          elements = "interfaces";
       };
    };

    disp = 
    {
      type = "LoadScale";

      scaleFunc = "0.0001 * (i)";

      model =
      {
        type     = "Constraints";
        conTable = "disp";
      };
    };

    lodi =
    {
      type  = "Lodi";
      group = "right";
      file  = "one-fiber-lodi.dat";
    };
  };
};
};

extraModules =
{
  modules = ["solver", "graph", "paraview", "view"];
  
  solver =  "Nonlin"
  {
      precision = 1.0e-5;
      maxIter   = 20;
  
      solver =
      {
        type = "SparseLU";
      };
  };

  paraview = "Paraview"
  {
    fileName = "one-fiber";
    elements = "bulk";
    dofs     = ["dx","dy"];
    data     = ["stress"];
  };

  graph =
  {
    type = "Graph";

    // Define the data sets to be visualized.
  
    dataSets = "loadDisp";
  
    loadDisp =
    {
      key = "Load-displacement curve";
      xData = "model.model.model.lodi.disp[0]";
      yData = "model.model.model.lodi.load[0]" ;
    };
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

    deformation = "50 * disp";

    // Define extra plugins for showing more data sets.

    plugins = "colors";

    colors =
    {
      type = "MeshColorView";
      data = "stress[stress_xx]";
      palette   = "custom";
      autoScale = true;
    };
  };

};
};

