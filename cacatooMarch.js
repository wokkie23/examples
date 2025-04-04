  if (typeof window == "undefined") {
    Simulation = require("../dist/cacatoo.js")
    }
  // NOTE: For a full version of this code with explanatory comments, see 0.html in the examples directory.
  // 1. SETUP. First, set up a configuration-object. Here we define how large the grid is, how long will it run, what colours will the critters be, etc.

  let sim
  var BirthRate = 0.9
  var DeathRate = 0.02
  var EnzymeProduction = 100
  var EnzymeDiffusion = 2
  var Mutation = 0.01
  var MutationExpression = 0.01
  var EnzymeDecayRate = 0.3
  var FoodLeft = 0.03
  var FoodRight = 0.03
  var Seperator = 1
  var GenesCost = 0.01

  var ExpressionCost = 0.1
  var expressionvisual = 0
  let logbook = ""
  let specTopLeft

  function cacatoo() {
  let config = {
    title: "cellular automatata 3 layers #1",
    description: "testing 2 squares",
    maxtime: 50000,
    ncol: 120,
    nrow: 120,
    scale: 2,
    graph_interval: 2000,  // Increased from 10
    graph_update: 2000,   // Increased from 100
    statecolours: {
      species: { 1: "darkgreen" },
      food: { 
        0: "black", 
        1: "#87CEFB",  // Short/processed
        2: "#4467C4",  // Medium
        3: "#00008C",  // Long
        4: "#FF6347",  // extra long
      },
      enzyme1: { 
        0: "black", 
        1: "#FFFF00",  // yellow
        2: "#FFD700",  // gold
        3: "#FFA500",  // orange
        4: "#FF8C00",  // dark orange
        5: "#FF4500",   // orange red
        6: "#FF0000",  // red
        7: "#DC143C",  // crimson
        8: "#8B0000",  // dark red
        9: "#400000",   // blood dark
        10: "#FF6347",  // tomato
        11: "#FF7F50",  // coral
        12: "#FF1493",  // deep pink
        13: "#FF69B4",  // hot pink
        14: "#FFB6C1",  // light pink
        15: "#FF00FF",  // magenta
        
      },
      enzyme2: {
        0: "#000000",
        1: "#FF7D7D",  // salmon pink
        2: "#FF4D4D",  // coral red
        3: "#D62222",  // fire engine red
        4: "#8B0000",  // dark red
        5: "#400000",   // blood dark
        6: "#FF0000",  // red
        7: "#DC143C",  // crimson
      },
      pget: {
        0: "black",
        1: "white",
        2: "red",
        3: "#0067BC",
        4: "#FF00FF",
        5: "grey",
        6: "#660000",
        7: "#000080",
        8: "purple",
        9: "#6A5ACD",
        10: "#35063E",
      },
      expression: { 0: "black", 1: "grey", 2: "darkgreen" },
      genes: { 0: "green", 1: "white", 2: "red", 3: "blue", 4: "#FF00FF" },
    },
    fpsmeter: true,
    wrap: [false, false],
    fastmode: false,
  }

  // 1. SETUP. (continued) Now, let's use that configuration-object to generate a new Cacatoo simulation
  sim = new Simulation(config)

  sim.makeGridmodel("master")

  let initialspecies = [
    { species: 1, genes: [0, 0], expression: [1, 1] },
    { species: 2 },
  ]

  let initialfood = [
    { food: { bitstring: "00000", value: 0 } },
    { food: { bitstring: "10001", value: 1 } },
    { food: { bitstring: "11011", value: 2 } },
    { food: { bitstring: "11111", value: 3 } }
  ]
  let enzyme1 = [{ enzyme1: 0 }]
  let enzyme2 = [{ enzyme2: 0 }]



  sim.populateGrid = function(gridmodel,individuals,freqs)
      {
          if(typeof gridmodel === 'string' || gridmodel instanceof String) gridmodel = this[gridmodel];
          if(individuals.length != freqs.length) throw new Error("populateGrid should have as many individuals as frequencies")
          if(freqs.reduce((a, b) => a + b) > 1) throw new Error("populateGrid should not have frequencies that sum up to greater than 1")

          for (let x = 0; x < gridmodel.nc; x++)                          // x are columns
              for (let y = 0; y < gridmodel.nr; y++){                 // y are rows
                  for (const property in individuals[0]) {
                      gridmodel.grid[x][y][property] = 0;    
                  }
                  let random_number = this.rng.random();
                  let sum_freqs = 0;
                  for(let n=0; n<individuals.length; n++)
                  {
                      sum_freqs += freqs[n];
                      if(random_number < sum_freqs) {
                          gridmodel.grid[x][y] = structuredClone(individuals[n]);
                          break
                      }
                  }
              }  
      }
      
  sim.populateGrid(sim.master, initialspecies, [0.5, 0])

  sim.makeGridmodel("particles")
  sim.makeGridmodel("enzyme1")
  sim.makeGridmodel("enzyme2")
  sim.createDisplay("master", "pget")
  //sim.createDisplay("master", "genes")
  sim.createDisplay("particles", "food")
  sim.createDisplay("enzyme1", "enzyme1")
  sim.createDisplay("enzyme2", "enzyme2")

  for (let x = 0; x < sim.enzyme1.nc; x++) {
    for (let y = 0; y < sim.enzyme1.nr; y++) {
      sim.enzyme1.grid[x][y] = { enzyme1: 0 };
      sim.enzyme2.grid[x][y] = { enzyme2: 0 };
    }
  }

  // 2. DEFINING THE RULES. Below, the user defines the nextState function. This function will be applied for each grid point when we will update the grid later.
  sim.master.nextState = function (i, j) {
    let RN = this.randomMoore8(this, i, j)
    let neighbours = sim.master.countMoore8(this, i, j, 1, "species")
    let state = this.grid[i][j].species
    let GeneA = this.grid[i][j].genes[0]
    let GeneB = this.grid[i][j].genes[1]
    let ExpressionA = this.grid[i][j].expression[0]
    let ExpressionB = this.grid[i][j].expression[1]

    if (this.grid[i][j].species == 0) {
      this.grid[i][j] = {
        species: 0,
        genes: [0, 0],
        expression: [0, 0],
      }
    }

    

    if (expressionvisual == true) {
      if (state === 0 || state === undefined) {
        this.grid[i][j].pget = 0;
      } else if (GeneA === 1 && GeneB === 1) {
        if (ExpressionA === 1 && ExpressionB === 1) {
          this.grid[i][j].pget = 4;
        } else if (ExpressionA === 1 && ExpressionB === 0) {
          this.grid[i][j].pget = 8;
        } else if (ExpressionB === 1 && ExpressionA === 0) {
          this.grid[i][j].pget = 9;
        } else if (ExpressionB === 0 && ExpressionA === 0){
          this.grid[i][j].pget = 10;
        }
      } else if (GeneA === 1) {
        this.grid[i][j].pget = ExpressionA === 1 ? 2 : 6;
      } else if (GeneB === 1) {
        this.grid[i][j].pget = ExpressionB === 1 ? 3 : 7;
      } else {
        this.grid[i][j].pget = (ExpressionA === 0 && ExpressionB === 0) ? 1 : 5;
      }
    } else {
      if (state === 0 || state === undefined) {
        this.grid[i][j].pget = 0;
      } else if (GeneA === 1 && GeneB === 1) {
        this.grid[i][j].pget = 4
      } else if (GeneA === 1) {
        this.grid[i][j].pget = 2
      } else if (GeneB === 1) {
        this.grid[i][j].pget = 3
      } else {
        this.grid[i][j].pget = 1
      }
    }


    if (state === 1 && this.rng.genrand_real1() < (EnzymeProduction / 100)) {
      const MAX_ENZYME = 6; 
      
      if (GeneA === 1 && GeneB === 1) {
          
          if (ExpressionA === 1 && ExpressionB === 1) {
              
              if (Math.random() < 0.5) {
                  if ((sim.enzyme1.grid[i][j].enzyme1 || 0) < MAX_ENZYME) {
                      sim.enzyme1.grid[i][j].enzyme1 = (sim.enzyme1.grid[i][j].enzyme1 || 0) + 1;
                  }
              } else {
                  if ((sim.enzyme2.grid[i][j].enzyme2 || 0) < MAX_ENZYME) {
                      sim.enzyme2.grid[i][j].enzyme2 = (sim.enzyme2.grid[i][j].enzyme2 || 0) + 1;
                  }
              }
          } else if (ExpressionA === 1) {
              // Only gene A expressed
              if ((sim.enzyme1.grid[i][j].enzyme1 || 0) < MAX_ENZYME) {
                  sim.enzyme1.grid[i][j].enzyme1 = (sim.enzyme1.grid[i][j].enzyme1 || 0) + 1;
              }
          } else if (ExpressionB === 1) {
              // Only gene B expressed
              if ((sim.enzyme2.grid[i][j].enzyme2 || 0) < MAX_ENZYME) {
                  sim.enzyme2.grid[i][j].enzyme2 = (sim.enzyme2.grid[i][j].enzyme2 || 0) + 1;
              }
          }
      } else if (GeneA === 1 && ExpressionA === 1) {
          // Only gene A present and expressed
          if ((sim.enzyme1.grid[i][j].enzyme1 || 0) < MAX_ENZYME) {
              sim.enzyme1.grid[i][j].enzyme1 = (sim.enzyme1.grid[i][j].enzyme1 || 0) + 1;
          }
      } else if (GeneB === 1 && ExpressionB === 1) {
          // Only gene B present and expressed
          if ((sim.enzyme2.grid[i][j].enzyme2 || 0) < MAX_ENZYME) {
              sim.enzyme2.grid[i][j].enzyme2 = (sim.enzyme2.grid[i][j].enzyme2 || 0) + 1;
          }
      }
  }

    let currentGenesCost
    let currentExpressionCost
    let effectiveBirthRate

    if (RN.genes[0] == 1 && RN.genes[1] == 1) {
      currentGenesCost = GenesCost * 2
    } else if (RN.genes[0] == 1 || RN.genes[1] == 1) {
      currentGenesCost = GenesCost
    } else if (RN.genes[0] == 0 && RN.genes[1] == 0) {
      currentGenesCost = 0
    }

    if (RN.expression[0] == 1 && RN.expression[1] == 1) {
      currentExpressionCost = ExpressionCost * 2
    } else if (RN.expression[0] == 1 || RN.expression[1] == 1) {
      currentExpressionCost = ExpressionCost
    } else if (RN.expression[0] == 0 && RN.expression[1] == 0) {
      currentExpressionCost = 0
    }

    effectiveBirthRate = BirthRate - currentGenesCost  - currentExpressionCost 

    if (
      this.grid[i][j].species == 0 &&
      RN.species == 1 &&
      sim.particles.grid[i][j].food == 1 &&
      this.rng.genrand_real1() < effectiveBirthRate
    ) {
      sim.particles.grid[i][j].food = 0

      this.grid[i][j] = {
        species: 1,
        genes: RN.genes,
        expression: RN.expression,
      }

      if (this.rng.genrand_real1() < Mutation) {
        this.grid[i][j] = {
          ...this.grid[i][j],
          genes: [
            this.rng.genrand_real1() < 0.5 ? 1 : 0,
            this.grid[i][j].genes[1],
          ],
        }
      }

      if (this.rng.genrand_real1() < Mutation) {
        this.grid[i][j] = {
          ...this.grid[i][j],
          genes: [
            this.grid[i][j].genes[0],
            this.rng.genrand_real1() < 0.5 ? 1 : 0,
          ],
        }
      }

      if (this.rng.genrand_real1() < MutationExpression) {
        this.grid[i][j] = {
          ...this.grid[i][j],
          expression: [
            this.rng.genrand_real1() < 0.5 ? 1 : 0,
            ExpressionB,
          ],
        }
      }

      if (this.rng.genrand_real1() < MutationExpression) {
        this.grid[i][j] = {
          ...this.grid[i][j],
          expression: [
            ExpressionA,
            this.rng.genrand_real1() < 0.5 ? 1 : 0,
          ],
        }
      }
    } else if (state == 1 && this.rng.genrand_real1() < DeathRate) {
      this.grid[i][j] = {
        ...this.grid[i][j],
        species: 0,
        genes: [0, 0],
        expression: [0, 0],
      }
    } else {
      this.grid[i][j].species = this.grid[i][j].species
    }
  }
  ///////////////////////////////////////////

  sim.particles.nextState = function (i, j) {
    if (sim.particles.grid[i][j].food == 2 && sim.enzyme1.grid[i][j].enzyme1 > 0) {
      sim.particles.grid[i][j].food = 1
    }
    if (sim.particles.grid[i][j].food == 3 && sim.enzyme2.grid[i][j].enzyme2 > 0) {
      sim.particles.grid[i][j].food = 2
    }
  }


  sim.enzyme1.nextState = function (i, j) {
    if (sim.enzyme1.grid[i][j].enzyme1 > 0) {
      if (this.rng.genrand_real1() < EnzymeDecayRate) {
        sim.enzyme1.grid[i][j].enzyme1 -= 1
      }
    }
  }
  sim.enzyme2.nextState = function (i, j) {
    if (sim.enzyme2.grid[i][j].enzyme2 > 0) {
      if (this.rng.genrand_real1() < EnzymeDecayRate) {
        sim.enzyme2.grid[i][j].enzyme2 -= 1
      }
    }
  }

  let influxFood1 = 0
  let influxFood2 = 0
  let influxFood3 = 0
  let influxFood4 = 0

  let foodProduced1 = 0
  let foodProduced2 = 0
  let foodProduced3 = 0

  let foodProduced1TopRight = 0
  let foodProduced2TopRight = 0
  let foodProduced3TopRight = 0

  let foodProduced1BottomRight = 0
  let foodProduced2BottomRight = 0
  let foodProduced3BottomRight = 0

  sim.particles.update = function () {
    foodProduced1TopRight = 0
    foodProduced2TopRight = 0
    foodProduced3TopRight = 0

    foodProduced1BottomRight = 0
    foodProduced2BottomRight = 0
    foodProduced3BottomRight = 0

    this.synchronous()
    let a = 0
    if (Seperator == 1) {
      // Ensure grid cells are initialized properly
      for (let i = 0; i < sim.particles.nc; i++) {
        for (let j = 0; j < sim.particles.nr; j++) {
          if (this.grid[i][j].food == undefined) {
            this.grid[i][j].food = 0
          }
        }
      }

      const spawnFood = (iStart, iEnd, jStart, jEnd, chance, foodType) => {
        for (let i = iStart; i < iEnd; i++) {
          for (let j = jStart; j < jEnd; j++) {
            if (
              this.grid[i][j].food == 0 &&
              this.rng.genrand_real1() < chance
            ) {
              this.grid[i][j].food = foodType

              const isTopHalf = i < sim.master.nc / 2
              const isLeftHalf = j < sim.master.nr / 2

              if (foodType == 1) {
                foodProduced1++
                if (!isTopHalf && isLeftHalf) {
                  foodProduced1TopRight++
                }
              } else if (foodType == 2) {
                foodProduced2++
                if (!isTopHalf && isLeftHalf) foodProduced2TopRight++
              } else if (foodType == 3) {
                foodProduced3++
                if (!isTopHalf && isLeftHalf) foodProduced3TopRight++
              }
            }
          }
        }
      }

      const spawnMixedFood = (iStart, iEnd, jStart, jEnd, chance) => {
        for (let i = iStart; i < iEnd; i++) {
          for (let j = jStart; j < jEnd; j++) {
            if (
              this.grid[i][j].food === 0 &&
              this.rng.genrand_real1() < chance
            ) {
              const rand = this.rng.genrand_real1()
              const isTopHalf = i < sim.master.nc / 2
              const isLeftHalf = j < sim.master.nr / 2

              if (rand < 0.3333) {
                this.grid[i][j].food = 1
                foodProduced1++
                if (!isTopHalf && !isLeftHalf) {
                  foodProduced1BottomRight++
                }
              } else if (rand < 0.6666) {
                this.grid[i][j].food = 2
                foodProduced2++
                if (!isTopHalf && !isLeftHalf) {
                  foodProduced2BottomRight++
                }
              } else if (rand < 1) {
                this.grid[i][j].food = 3
                foodProduced3++
                if (!isTopHalf && !isLeftHalf) {
                  foodProduced3BottomRight++
                }
              }
            }
          }
        }
      }

      const midCol = sim.particles.nc / 2
      const midRow = sim.particles.nr / 2

      // Top-left quadrant
      spawnFood(0, midCol, 0, midRow, FoodLeft, 1)

      // Bottom-left quadrant
      spawnFood(0, midCol, midRow, sim.particles.nr, FoodLeft, 2)

      // Right half
      spawnFood(midCol, sim.particles.nc, 0, midRow, FoodRight, 3)

      // Bottom-right quadrant with mixed food
      spawnMixedFood(
        midCol,
        sim.particles.nc,
        midRow,
        sim.particles.nr,
        FoodRight,
      )
    } else {
      for (let i = 0; i < sim.particles.nc; i++) {
        for (let j = 0; j < sim.particles.nr; j++) {
          if (this.grid[i][j].food == undefined) {
            this.grid[i][j].food = 0
          }
        }
      }

      const spawnFood = (iStart, iEnd, jStart, jEnd, chance, foodType) => {
        for (let i = iStart; i < iEnd; i++) {
          for (let j = jStart; j < jEnd; j++) {
            if (
              this.grid[i][j].food == 0 &&
              this.rng.genrand_real1() < chance
            ) {
              this.grid[i][j].food = foodType
            }
          }
        }
      }

      const midCol = sim.particles.nc / 2

      // Left half
      spawnFood(0, midCol, 0, sim.particles.nr, FoodLeft, 1)

      // Right half
      spawnFood(midCol, sim.particles.nc, 0, sim.particles.nr, FoodRight, 3)
    }
  }

  // Optimized update function with combined counting loops
  sim.master.update = function () {

    this.synchronous()
    let counts = {
      spec: {topLeft: 0, topRight: 0, bottomLeft: 0, bottomRight: 0},
      food: {1: 0, 2: 0, 3: 0},
      pget: {1: 0, 2: 0, 3: 0, 4: 0},
      mutations: {topLeft: 0, topRight: 0, bottomLeft: 0, bottomRight: 0}
    };
    
    // Single pass through the grid
    for (let x = 0; x < sim.master.nc; x++) {
      for (let y = 0; y < sim.master.nr; y++) {
        const area = getArea(x, y);
        
        // Species counting
        if (sim.master.grid[x][y].species == 1) {
          counts.spec[area]++;
        }
        
        // Food counting
        const foodType = sim.particles.grid[x][y].food;
        if (foodType > 0) {
          counts.food[foodType]++;
        }
        
        // Pget counting
        const pget = sim.master.grid[x][y].pget;
        if (pget >= 1 && pget <= 4) {
          counts.pget[pget]++;
        }
        
        // Mutation counting
        if (pget == 1 || pget == 2 || pget == 3 || pget == 4) {
          counts.mutations[area]++;
        }
      }
    }
    
    this.updateGraphs(counts.spec.bottomRight, counts.spec.topRight);
    
    function getArea(x, y) {
      const isTop = y < sim.master.nr / 2;
      const isLeft = x < sim.master.nc / 2;
      return isLeft ? (isTop ? 'topLeft' : 'bottomLeft') : (isTop ? 'topRight' : 'bottomRight');
    }
  }
  /////

  // New diffusion method that allows enzymes to spread further
  sim.enzyme1.diffuseMolecules = function() {
    // Create a temporary grid to track movements
    const newGrid = new Array(this.nc);
    for (let i = 0; i < this.nc; i++) {
      newGrid[i] = new Array(this.nr);
      for (let j = 0; j < this.nr; j++) {
        newGrid[i][j] = { enzyme1: 0 };
      }
    }

    // Process each cell
    for (let x = 0; x < this.nc; x++) {
      for (let y = 0; y < this.nr; y++) {
        const currentEnzymes = this.grid[x][y].enzyme1 || 0;
        
        // For each enzyme molecule in this cell
        for (let e = 0; e < currentEnzymes; e++) {
          // Decide whether to move (50% chance to stay)
          if (this.rng.genrand_real1() < 0) {
            newGrid[x][y].enzyme1 = (newGrid[x][y].enzyme1 || 0) + 1;
            continue;
          }

          // Choose a random neighbor (Moore neighborhood)
          let nx = x, ny = y;
          const dir = Math.floor(this.rng.genrand_real1() * 8);
          switch(dir) {
            case 0: nx = x-1; ny = y-1; break;
            case 1: nx = x;   ny = y-1; break;
            case 2: nx = x+1; ny = y-1; break;
            case 3: nx = x-1; ny = y;   break;
            case 4: nx = x+1; ny = y;   break;
            case 5: nx = x-1; ny = y+1; break;
            case 6: nx = x;   ny = y+1; break;
            case 7: nx = x+1; ny = y+1; break;
          }

          // Wrap around or stay at edge based on simulation settings
          if (config.wrap[0]) {
            nx = (nx + this.nc) % this.nc;
          } else {
            nx = Math.max(0, Math.min(this.nc-1, nx));
          }
          
          if (config.wrap[1]) {
            ny = (ny + this.nr) % this.nr;
          } else {
            ny = Math.max(0, Math.min(this.nr-1, ny));
          }

          newGrid[nx][ny].enzyme1 = (newGrid[nx][ny].enzyme1 || 0) + 1;
        }
      }
    }

    // Copy back to main grid
    for (let i = 0; i < this.nc; i++) {
      for (let j = 0; j < this.nr; j++) {
        this.grid[i][j].enzyme1 = newGrid[i][j].enzyme1;
      }
    }
  };

  // Similar implementation for enzyme2
  sim.enzyme2.diffuseMolecules = function() {
    const newGrid = new Array(this.nc);
    for (let i = 0; i < this.nc; i++) {
      newGrid[i] = new Array(this.nr);
      for (let j = 0; j < this.nr; j++) {
        newGrid[i][j] = { enzyme2: 0 };
      }
    }

    for (let x = 0; x < this.nc; x++) {
      for (let y = 0; y < this.nr; y++) {
        const currentEnzymes = this.grid[x][y].enzyme2 || 0;
        
        for (let e = 0; e < currentEnzymes; e++) {
          if (this.rng.genrand_real1() < 0.5) {
            newGrid[x][y].enzyme2 = (newGrid[x][y].enzyme2 || 0) + 1;
            continue;
          }

          let nx = x, ny = y;
          const dir = Math.floor(this.rng.genrand_real1() * 8);
          switch(dir) {
            case 0: nx = x-1; ny = y-1; break;
            case 1: nx = x;   ny = y-1; break;
            case 2: nx = x+1; ny = y-1; break;
            case 3: nx = x-1; ny = y;   break;
            case 4: nx = x+1; ny = y;   break;
            case 5: nx = x-1; ny = y+1; break;
            case 6: nx = x;   ny = y+1; break;
            case 7: nx = x+1; ny = y+1; break;
          }

          if (config.wrap[0]) {
            nx = (nx + this.nc) % this.nc;
          } else {
            nx = Math.max(0, Math.min(this.nc-1, nx));
          }
          
          if (config.wrap[1]) {
            ny = (ny + this.nr) % this.nr;
          } else {
            ny = Math.max(0, Math.min(this.nr-1, ny));
          }

          newGrid[nx][ny].enzyme2 = (newGrid[nx][ny].enzyme2 || 0) + 1;
        }
      }
    }

    for (let i = 0; i < this.nc; i++) {
      for (let j = 0; j < this.nr; j++) {
        this.grid[i][j].enzyme2 = newGrid[i][j].enzyme2;
      }
    }
  };

  // Update the update functions to use the new diffusion method
  sim.enzyme1.update = function () {
    this.synchronous();
    if (this.time % EnzymeDiffusion == 0) {
      this.diffuseMolecules();
    }
    
    // Enzyme decay
    for (let i = 0; i < this.nc; i++) {
      for (let j = 0; j < this.nr; j++) {
        if (this.grid[i][j].enzyme1 > 0 && this.rng.genrand_real1() < EnzymeDecayRate) {
          this.grid[i][j].enzyme1 -= 1;
        }
      }
    }
  };

  sim.enzyme2.update = function () {
    this.synchronous();
    if (this.time % EnzymeDiffusion == 0) {
      this.diffuseMolecules();
    }
    
    // Enzyme decay
    for (let i = 0; i < this.nc; i++) {
      for (let j = 0; j < this.nr; j++) {
        if (this.grid[i][j].enzyme2 > 0 && this.rng.genrand_real1() < EnzymeDecayRate) {
          this.grid[i][j].enzyme2 -= 1;
        }
      }
    }
  };

  /////

  sim.addSlider("expressionvisual", 0, 1, 1)
  sim.addButton("step", function () {
    sim.step()
    sim.display()
  })

  sim.addButton("Well mix", function () {
    sim.toggle_mix()
  })

  sim.addButton("pauze", function () {
    sim.toggle_play()
  })

  sim.master.updateGraphs = function(specBottomRight, specTopRight) {
    if(this.time>=4000) {
      if (this.time % 100 == 0) {
        let log = `time: ${this.time}, bottom: ${specBottomRight}, top: ${specTopRight}`
        sim.log(log, "output")
        logbook += log + ''
      }
    }
    
    // if(this.time==15000) {
    //   sim.write(logbook,"bottomvstop.txt")
    // }
  }

  function mixQuadrant(gridmodel, quadrant) {
    function mixSingleGrid(grid) {
        const midCol = Math.floor(grid.nc / 2);
        const midRow = Math.floor(grid.nr / 2);
        let startX, endX, startY, endY;
        switch (quadrant) {
            case "top-left":
                startX = 0; endX = midCol; startY = 0; endY = midRow; break;
            case "top-right":
                startX = midCol; endX = grid.nc; startY = 0; endY = midRow; break;
            case "bottom-left":
                startX = 0; endX = midCol; startY = midRow; endY = grid.nr; break;
            case "bottom-right":
                startX = midCol; endX = grid.nc; startY = midRow; endY = grid.nr; break;
            default:
                throw new Error("Invalid quadrant specified");
        }
        let cells = [];
        for (let x = startX; x < endX; x++) {
            for (let y = startY; y < endY; y++) {
                if (Seperator === 1) {
                    cells.push(grid.grid[x][y]);
                }
            }
        }
        for (let i = cells.length - 1; i > 0; i--) {
            const j = Math.floor(sim.rng.genrand_real1() * (i + 1));
            [cells[i], cells[j]] = [cells[j], cells[i]];
        }
        let index = 0;
        for (let x = startX; x < endX; x++) {
            for (let y = startY; y < endY; y++) {
                if (Seperator === 1) {
                    grid.grid[x][y] = cells[index++];
                }
            }
        }
    }

    mixSingleGrid(gridmodel);
    mixSingleGrid(sim.particles);
    mixSingleGrid(sim.enzyme1);
    mixSingleGrid(sim.enzyme2);
  }

  sim.addButton("Mix Bottom Right", function () {
    mixBottomRightQuadrant(sim.master);
    sim.display();
  });

  sim.addSlider("Seperator", 0, 1, 1)
  sim.addSlider("BirthRate", 0, 1, 0.1)
  sim.addSlider("EnzymeDecayRate", 0, 1, 0.01)
  sim.addSlider("EnzymeProduction", 0, 100, 1)
  sim.addSlider("DeathRate", 0, 0.1, 0.005)
  sim.addSlider("EnzymeDiffusion", 1, 50, 1)
  sim.addSlider("Mutation", 0, 1, 0.01)
  sim.addSlider("FoodLeft", 0, 0.1, 0.01)
  sim.addSlider("FoodRight", 0, 0.1, 0.01)
  sim.addSlider("ExpressionCost", 0, 1, 0.01)
  sim.addSlider("GenesCost", 0, 1, 0.05)

  sim.start()
  sim.addStatebrush(sim.master, "species", 0, 200, 1)
  }

  if (typeof window == "undefined") cacatoo()