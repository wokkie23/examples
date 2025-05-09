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
    var EnzymeDecayRate = 0.35
    var FoodDecayRate = 0
    var FoodLeft = 0.02
    var FoodRight = 0.02
    var Seperator = 1
    var GenesCost = 0.01

    var ExpressionCost = 0.1
    var expressionvisual = 0
    let logbook = ""
    let specTopLeft

    function cacatoo() {
    let config = {
      title: "cellular automatata 3 layers #1",
      description: "",
      maxtime: 50000,
      ncol: 120,
      nrow: 120,
      scale: 2,
      graph_interval: 20000,  
      graph_update: 20000,   
      statecolours: {
        species: { 1: "darkgreen" },
        displayValue: {  
          0: "black",        // 0: No food
          1: "#FF0000",     // 1: Vivid Red (length 1)
          2: "#DC143C",     // 2: Crimson (length 2)
          3: "#FF00FF",     // 3: Magenta (length 3)
          4: "purple",     // 4: Orange Red (length 4)
          5: "#FFA500",     // 5: Orange (length 5)
          6: "#FFD700",     // 6: Gold (length 6)
          7: "#FFFF00",     // 7: Yellow (length 7)
          8: "#ADFF2F",     // 8: Green-Yellow (length 8)
          9: "#7FFF00",     // 9: Chartreuse (length 9)
          10: "#00FF00",    // 10: Pure Green (length 10)
          11: "#00FA9A",    // 11: Medium Spring Green
          12: "#00CED1",    // 12: Dark Turquoise
          13: "#1E90FF",    // 13: Dodger Blue
          14: "#0000FF",    // 14: Pure Blue
          15: "#8A2BE2"     // 15: Blue Violet (longest)
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
          15: "#FF00FF",
          16: "#FF00FF",
          17: "#FF00FF",
          18: "#FF00FF",
          19: "#FF00FF",  // magenta
          20: "#FF00FF",
          21: "#FF00FF",
          22: "#FF00FF",
          23: "#FF00FF",
          24: "#FF00FF"
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
          8: "#DC143C",
          9: "#DC143C"
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
      },
      fpsmeter: true, 
      wrap: [false, false],
      fastmode: true,
    }

    // 1. SETUP. (continued) Now, let's use that configuration-object to generate a new Cacatoo simulation
    sim = new Simulation(config)

    sim.makeGridmodel("master")

    let initialspecies = [
      { species: 1, genes: [0, 0], expression: [1, 1] },
      { species: 2 },
    ]

    var FoodConfig = {
      upleft: ["1", "1"],    // Default left patterns
      upright: ["111", "111"], // Default right patterns
      downleft: ["101", "10101"],
      downright: ["101011111", "101011111"]
    };

    var MaxBitLength = 10



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

    for (let x = 0; x < sim.particles.nc; x++) {
      for (let y = 0; y < sim.particles.nr; y++) {
          sim.particles.grid[x][y] = { 
              food: [],  // Start with empty array
              displayValue: 0 
          };
      }
    }

    sim.makeGridmodel("enzyme1")
    sim.makeGridmodel("enzyme2")
    sim.createDisplay("master", "pget")
    //sim.createDisplay("master", "genes")
    // sim.createDisplay("particles", "food")
    sim.createDisplay("particles", "displayValue")
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

      if (sim.particles.grid[i][j].food.bitstring) {
        const ones = (sim.particles.grid[i][j].food.bitstring.match(/1/g) || []).length;
        const total = sim.particles.grid[i][j].food.bitstring.length || 1;
        const ratio = ones / total;
        sim.particles.grid[i][j].food.value = Math.floor(ratio * 15) || 0;
      }


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


    //   if (state === 1 && this.rng.genrand_real1() < (EnzymeProduction / 100)) {
    //     const MAX_ENZYME = 6; 
        
    //     if (GeneA === 1 && GeneB === 1) {
            
    //         if (ExpressionA === 1 && ExpressionB === 1) {
                
    //             if (Math.random() < 0.5) {
    //                 if ((sim.enzyme1.grid[i][j].enzyme1 || 0) < MAX_ENZYME) {
    //                     sim.enzyme1.grid[i][j].enzyme1 = (sim.enzyme1.grid[i][j].enzyme1 || 0) + 1;
    //                 }
    //             } else {
    //                 if ((sim.enzyme2.grid[i][j].enzyme2 || 0) < MAX_ENZYME) {
    //                     sim.enzyme2.grid[i][j].enzyme2 = (sim.enzyme2.grid[i][j].enzyme2 || 0) + 1;
    //                 }
    //             }
    //         } else if (ExpressionA === 1) {
    //             // Only gene A expressed
    //             if ((sim.enzyme1.grid[i][j].enzyme1 || 0) < MAX_ENZYME) {
    //                 sim.enzyme1.grid[i][j].enzyme1 = (sim.enzyme1.grid[i][j].enzyme1 || 0) + 1;
    //             }
    //         } else if (ExpressionB === 1) {
    //             // Only gene B expressed
    //             if ((sim.enzyme2.grid[i][j].enzyme2 || 0) < MAX_ENZYME) {
    //                 sim.enzyme2.grid[i][j].enzyme2 = (sim.enzyme2.grid[i][j].enzyme2 || 0) + 1;
    //             }
    //         }
    //     } else if (GeneA === 1 && ExpressionA === 1) {
    //         // Only gene A present and expressed
    //         if ((sim.enzyme1.grid[i][j].enzyme1 || 0) < MAX_ENZYME) {
    //             sim.enzyme1.grid[i][j].enzyme1 = (sim.enzyme1.grid[i][j].enzyme1 || 0) + 1;
    //         }
    //     } else if (GeneB === 1 && ExpressionB === 1) {
    //         // Only gene B present and expressed
    //         if ((sim.enzyme2.grid[i][j].enzyme2 || 0) < MAX_ENZYME) {
    //             sim.enzyme2.grid[i][j].enzyme2 = (sim.enzyme2.grid[i][j].enzyme2 || 0) + 1;
    //         }
    //     }
    // }

    if (state === 1 && this.rng.genrand_real1() < (EnzymeProduction / 100)) {
      // Enzyme1 produces "11" pattern
      if (GeneA === 1 && ExpressionA === 1) {
        sim.enzyme1.grid[i][j].enzyme1 = (sim.enzyme1.grid[i][j].enzyme1 || 0) + 1;
      }
      // Enzyme2 could produce a different pattern
      if (GeneB === 1 && ExpressionB === 1) {
        sim.enzyme2.grid[i][j].enzyme2 = (sim.enzyme2.grid[i][j].enzyme2 || 0) + 1;
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

      function hasLengthOneFood(foods) {
        return foods.some(food => food.bitstring.length === 1);
    }

      if (
        this.grid[i][j].species == 0 &&
        RN.species == 1 &&
        sim.particles.grid[i][j].food.some(food => food.bitstring === "1") && // Check if any food exists
        this.rng.genrand_real1() < effectiveBirthRate
    ) {
      // Remove FIRST length-1 food item found
      const foods = sim.particles.grid[i][j].food;
      const foodIndex = foods.findIndex(food => food.bitstring.length === 1);
      if(foodIndex > -1) {
        sim.particles.grid[i][j].food.splice(foodIndex, 1);
    }


        //sim.particles.grid[i][j].food = 0 

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
      // Convert food storage to array format if needed
      if (!Array.isArray(this.grid[i][j].food)) {
          this.grid[i][j].food = this.grid[i][j].food.bitstring ? [this.grid[i][j].food] : [];
      }

      this.grid[i][j].food = this.grid[i][j].food.filter(foodItem => {
        // Explicit version checking both length and content
        return !(foodItem.bitstring.length === 1 && foodItem.bitstring === "0");
      });



      // Process each food item in the cell
      const newFoods = [];
      let didSplit = false;

      for (const foodItem of this.grid[i][j].food) {
          if (!foodItem.bitstring) continue;

          let processed = false;
          
          // Check enzyme1 action (split "11" patterns)
          if (sim.enzyme1.grid[i][j].enzyme1 > 0) {
            const pattern = "10";
            const index = foodItem.bitstring.indexOf(pattern);
        
            if (index === 0) { // Only split if "11" is at the start
                const left = foodItem.bitstring.substring(0, pattern.length);
                const right = foodItem.bitstring.substring(pattern.length);
        
                // Split "11" into two "1"s
                for (let i = 0; i < left.length; i++) {
                    newFoods.push({
                        bitstring: left[i],
                        value: 0.5 // Since "1" is 1/2 in value
                    });
                }
        
                // If there's remaining string, add it
                if (right.length > 0) {
                    newFoods.push({
                        bitstring: right,
                        value: parseInt(right, 2) / (1 << right.length)
                    });
                }
        
                // IMPORTANT: Mark the original food item for removal
                processed = true;
                didSplit = true;
            }
        }

        if (sim.enzyme2.grid[i][j].enzyme2 > 0) {
          const pattern = "11";
          const index = foodItem.bitstring.indexOf(pattern);
      
          if (index === 0) { // Only split if "11" is at the start
              const left = foodItem.bitstring.substring(0, pattern.length);
              const right = foodItem.bitstring.substring(pattern.length);
      
              // Split "11" into two "1"s
              for (let i = 0; i < left.length; i++) {
                  newFoods.push({
                      bitstring: left[i],
                      value: 0.5 // Since "1" is 1/2 in value
                  });
              }
      
              // If there's remaining string, add it
              if (right.length > 0) {
                  newFoods.push({
                      bitstring: right,
                      value: parseInt(right, 2) / (1 << right.length)
                  });
              }
      
              // IMPORTANT: Mark the original food item for removal
              processed = true;
              didSplit = true;
          }
      }

          if (!processed) {
              // Keep original food if not processed
              newFoods.push(foodItem);
          }
      }

      // Update food array for the cell
      this.grid[i][j].food = newFoods;

      let maxLength = 0;
      this.grid[i][j].food.forEach(food => {
          const len = food.bitstring.length;
          if (len > maxLength) maxLength = len;
      });
      
      // Scale length to 0-15 range based on MaxBitLength
      const lengthValue = Math.min(15, Math.floor((maxLength / MaxBitLength) * 15));
      this.grid[i][j].displayValue = lengthValue;
        
    };


    sim.enzyme1.nextState = function(i, j) {
      let enzyme1 = sim.enzyme1.grid[i][j].enzyme1;
      if (enzyme1 > 0) {
        // Probabilistically remove each enzyme (30% chance per unit)
        let removed = 0;
        for (let i = 0; i < enzyme1; i++) {
          if (this.rng.genrand_real1() < EnzymeDecayRate) removed++;
        }
        sim.enzyme1.grid[i][j].enzyme1 = enzyme1 - removed;
      }
    };
    sim.enzyme2.nextState = function(i, j) {
      let enzyme2 = sim.enzyme2.grid[i][j].enzyme2;
      if (enzyme2 > 0) {
        // Probabilistically remove each enzyme (30% chance per unit)
        let removed = 0;
        for (let i = 0; i < enzyme2; i++) {
          if (this.rng.genrand_real1() < EnzymeDecayRate) removed++;
        }
        sim.enzyme2.grid[i][j].enzyme2 = enzyme2 - removed;
      }
    };

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

      for (let i = 0; i < sim.particles.nc; i++) {
        for (let j = 0; j < sim.particles.nr; j++) {
            if (!sim.particles.grid[i][j].food) {
                sim.particles.grid[i][j].food = { bitstring: "", value: 0 };
            }
        }
    }

      this.asynchronous()
      let a = 0
      if (Seperator == 1) {
        // Ensure grid cells are initialized properly
        for (let i = 0; i < sim.particles.nc; i++) {
          for (let j = 0; j < sim.particles.nr; j++) {
            if (this.grid[i][j].food == undefined) {
              //this.grid[i][j].food = 0
              this.grid[i][j].food = { bitstring: "", value: 0 };
            }
          }
        }

        // const spawnFood = (iStart, iEnd, jStart, jEnd, chance, foodType) => {
        //   for (let i = iStart; i < iEnd; i++) {
        //     for (let j = jStart; j < jEnd; j++) {
        //       if (
        //         this.grid[i][j].food == 0 &&
        //         this.rng.genrand_real1() < chance
        //       ) {
        //         this.grid[i][j].food = foodType

        //         const isTopHalf = i < sim.master.nc / 2
        //         const isLeftHalf = j < sim.master.nr / 2

        //         if (foodType == 1) {
        //           foodProduced1++
        //           if (!isTopHalf && isLeftHalf) {
        //             foodProduced1TopRight++
        //           }
        //         } else if (foodType == 2) {
        //           foodProduced2++
        //           if (!isTopHalf && isLeftHalf) foodProduced2TopRight++
        //         } else if (foodType == 3) {
        //           foodProduced3++
        //           if (!isTopHalf && isLeftHalf) foodProduced3TopRight++
        //         }
        //       }
        //     }
        //   }
        // }

        // const spawnMixedFood = (iStart, iEnd, jStart, jEnd, chance) => {
        //   for (let i = iStart; i < iEnd; i++) {
        //     for (let j = jStart; j < jEnd; j++) {
        //       if (
        //         this.grid[i][j].food === 0 &&
        //         this.rng.genrand_real1() < chance
        //       ) {
        //         const rand = this.rng.genrand_real1()
        //         const isTopHalf = i < sim.master.nc / 2
        //         const isLeftHalf = j < sim.master.nr / 2

        //         if (rand < 0.3333) {
        //           this.grid[i][j].food = 1
        //           foodProduced1++
        //           if (!isTopHalf && !isLeftHalf) {
        //             foodProduced1BottomRight++
        //           }
        //         } else if (rand < 0.6666) {
        //           this.grid[i][j].food = 2
        //           foodProduced2++
        //           if (!isTopHalf && !isLeftHalf) {
        //             foodProduced2BottomRight++
        //           }
        //         } else if (rand < 1) {
        //           this.grid[i][j].food = 3
        //           foodProduced3++
        //           if (!isTopHalf && !isLeftHalf) {
        //             foodProduced3BottomRight++
        //           }
        //         }
        //       }
        //     }
        //   }
        // }

        const midCol = sim.particles.nc / 2
        const midRow = sim.particles.nr / 2

        const spawnBitstringFood = (iStart, iEnd, jStart, jEnd, chance, patterns) => {
          for (let i = iStart; i < iEnd; i++) {
            for (let j = jStart; j < jEnd; j++) {
              if ((sim.particles.grid[i][j].food.length < 1 && // Limit to 2 food items
                this.rng.genrand_real1() < chance)) {
                // Select random pattern from available ones
                const pattern = patterns[Math.floor(this.rng.genrand_real1() * patterns.length)];
                const ones = (pattern.match(/1/g) || []).length;
                
                sim.particles.grid[i][j].food.push({
                  bitstring: pattern,
                  value: ones / pattern.length
                });
              }
            }
          }
        };
        
        spawnBitstringFood(0, midCol, 0, midRow, FoodLeft, FoodConfig.upleft);
        spawnBitstringFood(0, midCol, midRow, sim.particles.nr, FoodLeft, FoodConfig.downleft);
        spawnBitstringFood(midCol, sim.particles.nc, 0, midRow, FoodRight, FoodConfig.upright);
        spawnBitstringFood(midCol, sim.particles.nc, midRow, sim.particles.nr, FoodRight, FoodConfig.downright);

        for (let i = 0; i < this.nc; i++) {
          for (let j = 0; j < this.nr; j++) {
            if (sim.particles.grid[i][j].food.length > 0 && this.rng.genrand_real1() < FoodDecayRate) {
              sim.particles.grid[i][j].food = [];
              sim.particles.grid[i][j].displayValue = 0
            }
          }
        }
    } }

    // Optimized update function with combined counting loops
    sim.master.update = function () {

      this.asynchronous()
      let counts = {
        spec: {topLeft: 0, topRight: 0, bottomLeft: 0, bottomRight: 0},
        food: {1: 0, 2: 0, 3: 0},
        pget: {1: 0, 2: 0, 3: 0, 4: 0},
        mutations: {topLeft: 0, topRight: 0, bottomLeft: 0, bottomRight: 0}
      };
      
      // // Single pass through the grid
      // for (let x = 0; x < sim.master.nc; x++) {
      //   for (let y = 0; y < sim.master.nr; y++) {
      //     const area = getArea(x, y);
          
      //     // Species counting
      //     if (sim.master.grid[x][y].species == 1) {
      //       counts.spec[area]++;
      //     }
          
      //     // Food counting
      //     const foodType = sim.particles.grid[x][y].food;
      //     if (foodType > 0) {
      //       counts.food[foodType]++;
      //     }
          
      //     // Pget counting
      //     const pget = sim.master.grid[x][y].pget;
      //     if (pget >= 1 && pget <= 4) {
      //       counts.pget[pget]++;
      //     }
          
      //     // Mutation counting
      //     if (pget == 1 || pget == 2 || pget == 3 || pget == 4) {
      //       counts.mutations[area]++;
      //     }
      //   }
      // }
      
      // this.updateGraphs(counts.spec.bottomRight, counts.spec.topRight);
      
      // function getArea(x, y) {
      //   const isTop = y < sim.master.nr / 2;
      //   const isLeft = x < sim.master.nc / 2;
      //   return isLeft ? (isTop ? 'topLeft' : 'bottomLeft') : (isTop ? 'topRight' : 'bottomRight');
      // }
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
      this.asynchronous();
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
      this.asynchronous();
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
    sim.addSlider("MutationExpression", 0, 1, 0.01)

    sim.start()
    sim.addStatebrush(sim.master, "species", 0, 200, 1)
    }

    if (typeof window == "undefined") cacatoo()
      