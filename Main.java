package FrogsLinear2024_3;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;

import TwoDimensions_Monopolization.Individual;

public class Main {

	public static void main(String[] args) throws IOException {
		
		int maxGen = 1000;
		int initialPopSize = 2000;
		int chromosomeLength = 100;
		int matingRadius = 1;
		double similarityThreshold = 0.95;
		double v = 0.001; 
		double m = 0.001;

		Random rand = new Random();
		Space map = new Space(chromosomeLength);
		
		
		File file = new File("/Users/feper/eclipse-workspace/PolyploidyDeterministicEqs/Monopolization8.txt");
		FileWriter fw = new FileWriter(file);
		PrintWriter speciesFile = new PrintWriter(fw);
		speciesFile.println("t;clustering;niche_4x;niche_2x");
		
		int sample = 0;
		while(sample < 1) {
			
			ArrayList<Individual> population = new ArrayList<>();
			
			while(population.size() < initialPopSize) {
				
				int[][] genotype = new int[2][chromosomeLength];
				
				int posX = (int) (Math.random() * (100));
		        int posY = (int) (Math.random() * (100));
		        
		        int key = map.getMap()[0][0];
		        if(key > 0) {

		        	int[] reference = map.getRefGenome(key);
			        
					for(int j = 0; j < reference.length; j++) {
						genotype[0][j] = reference[j];
						genotype[1][j] = reference[j];
					}
					population.add(new Individual(matingRadius, posX, posY, m, v, genotype, reference));
					
		        }
			}
			
			int t = 1;
			while(t <= maxGen) {
				
				ArrayList<Individual> newPop = new ArrayList<>();
				Collections.shuffle(population);
				
				for(int i = 0; i < population.size(); i++) {
					
					int index = i;
					Individual mate01 = population.get(index);
					
					ArrayList<Individual> neighbours = scanNeighbourhood(mate01, index, population, matingRadius, map.getMap());
					if(neighbours.size() > 0) {
						
						Individual mate02 = neighbours.get(rand.nextInt(neighbours.size()));
						
						int n_offspring = getPoissonRandom(9);
						if(mate01.getGenotype().length == 4) {
							n_offspring = getPoissonRandom(4);
						}
						//if(n_offspring < 0) n_offspring = 1;
						for(int j = 0; j < n_offspring; j++) {
							
							try {
								
								int[][] g1 = mate01.getGametes();
								int[][] g2 = mate02.getGametes();
								
								boolean samePloidy = true;
								/*if(g1.length == g2.length) {
									samePloidy = true;
								}*/
								
								if(samePloidy) {
									
									//double similarity = compare(g1[0], g2[0]);
									//if(similarity > similarityThreshold) {
										
										int[] positions = newPositions(mate01.getX_coord(), mate01.getY_coord(), matingRadius, map.getMap());
										int offspring_x = positions[0];
										int offspring_y = positions[1];
										
										int[][] offspring = syngamy(g1, g2);
										int offspringKey = map.getMap()[offspring_x][offspring_y];
										int[] offspringReference = map.getRefGenome(offspringKey);
										
										double m1 = m;
										/*if(offspring.length == 4) {
											m1 = 10.0*m;
										}*/
										Individual ind = new Individual(matingRadius, offspring_x, offspring_y, m1, v, offspring, offspringReference);
										ArrayList<Individual> neighboursOffspring = scanNeighbourhood(ind, -1, newPop, matingRadius, map.getMap());
										double n = (double) neighboursOffspring.size();
										double cells = 9.0;
										double density = n/cells;
										
										double f1 = computeFitness(ind.getGenotype(), offspringReference);
										//double select = Math.random();
										//System.out.println(f1);
										if(density <= f1) {
											
											newPop.add(ind);				
										}else {
											break;
										}
										
									//}
								}

							}catch(NullPointerException e) {
								
							}
						}
					}
				}
				
				population = newPop;
				
				File file2 = new File("/Users/feper/eclipse-workspace/PolyploidyDeterministicEqs/Monopolization8/Time"+t+".txt");
				FileWriter fw2 = new FileWriter(file2);
				PrintWriter polyploidMap2 = new PrintWriter(fw2);
				polyploidMap2.println("x_coord;y_coord;ploidy;fitness");
				
				for(int i = 0; i < population.size(); i++) {
					
					int key = map.getMap()[population.get(i).getX_coord()][population.get(i).getY_coord()];
					int[] keyRef = map.getRefGenome(key);
					double fitness = computeFitness(population.get(i).getGenotype(), keyRef);
					
					polyploidMap2.println(population.get(i).getX_coord()+";"+population.get(i).getY_coord()+";"+population.get(i).getPloidy()+";"+fitness);
				}
				
				polyploidMap2.close();

				double generalFrequency = 0;
				double countTetraploids = 0;
				double countTriploids = 0;
				double countDiploids = 0;
				double nicheTetraploids = 0;
				double nicheTriploids = 0;
				double nicheDiploids = 0;
				
				for(int i = 0; i < population.size(); i++) {
					Individual ind = population.get(i);
					if(ind.getPloidy() == 4) {
						
						countTetraploids += 1.0;
						nicheTetraploids += map.getMap()[ind.getX_coord()][ind.getY_coord()];
						
						ArrayList<Individual> neighbours = scanNeighbourhood(ind, i, population, matingRadius, map.getMap());
						
						double numberTetra = 0.0;
						
						for(int j = 0; j < neighbours.size(); j++) {
							Individual test = neighbours.get(j);
							if(test.getPloidy() == 4) {
								numberTetra += 1.0;
							}
						}
						
						if(neighbours.size() > 0 ) {
							generalFrequency += numberTetra/neighbours.size(); 
						}
							
					}else if(ind.getPloidy() == 2) {
						
						countDiploids += 1.0;
						nicheDiploids += map.getMap()[ind.getX_coord()][ind.getY_coord()];
						
					}
				}
				
				if(countTetraploids > 0) {
					generalFrequency = generalFrequency/countTetraploids;
					nicheTetraploids = nicheTetraploids/countTetraploids;
				}else if(countDiploids > 0){
					generalFrequency = 0.0;
					nicheDiploids = nicheDiploids/countDiploids;
				}
					

				if(countDiploids > 0) {
					nicheDiploids = nicheDiploids/countDiploids;
				}
				
				System.out.println("Population size: "+population.size()+"\t Generation: "+t+"\t"+countTetraploids/population.size()+"\t"+generalFrequency);
				
				speciesFile.println(t+";"+generalFrequency+";"+nicheTetraploids+";"+nicheDiploids);

				t++;
			}
			sample++;
		}
		
		
		
		speciesFile.close();
		
	}

	public static double compare(int[] g1, int[] g2) {
		
		double c = 0.0;
		for(int i = 0; i < g1.length; i++) {
			if(g1[i] == g2[i]) {
				c += 1.0;
			}
		}
		
		return c/g1.length;
	}
	
	
	public static ArrayList<Individual> scanNeighbourhood(Individual mate01, int index, ArrayList<Individual> population, int radius, int[][] map) {
		
		ArrayList<Individual> neighbours = new ArrayList<>();
		
		int x = mate01.getX_coord();
		int y = mate01.getY_coord();
		int start_x = x - radius, end_x = x + radius;
		int start_y = y - radius, end_y = y + radius;
		
		if(start_x < 0) {
			start_x = 199 + start_x;
		}
		
		if(end_x > 199) {
			end_x = 0 + (end_x - 199);
		}
		
		if(start_y < 0) {
			start_y = 99 + start_y;
		}
		
		if(end_y > 99) {
			end_y = 0 + (end_y - 99);
		}
		
		for(int k = 0; k < population.size(); k++) {
			
			if(start_x < end_x && start_y < end_y) {
				
				if(population.get(k).getX_coord() >= start_x && population.get(k).getX_coord() <= end_x
						&& population.get(k).getY_coord() >= start_y && population.get(k).getY_coord() <= end_y &&
						k != index) {
					
					neighbours.add(population.get(k));
					
				}
				
			}else if(start_x > end_x && start_y < end_y) {
				
				if(population.get(k).getX_coord() >= start_x && population.get(k).getX_coord() <= 199 
						|| population.get(k).getX_coord() >= 0 && population.get(k).getX_coord() <= end_x
						&& population.get(k).getY_coord() >= start_y && population.get(k).getY_coord() <= end_y &&
						k != index) {
					
					neighbours.add(population.get(k));
					
				}
				
			}else if(start_x < end_x && start_y > end_y) {
				
				if(population.get(k).getX_coord() >= start_x && population.get(k).getX_coord() <= end_x
						&& population.get(k).getY_coord() >= start_y && population.get(k).getY_coord() <= 99 
						|| population.get(k).getY_coord() >= 0 && population.get(k).getY_coord() <= end_y 
						&& k != index) {
					
					neighbours.add(population.get(k));
					
				}
				
			}else if(start_x > end_x && start_y > end_y) {
				
				if(population.get(k).getX_coord() >= start_x && population.get(k).getX_coord() <= 199 
						|| population.get(k).getX_coord() >= 0 && population.get(k).getX_coord() <= end_x
						&& population.get(k).getY_coord() >= start_y && population.get(k).getY_coord() <= 99 
						|| population.get(k).getY_coord() >= 0 && population.get(k).getY_coord() <= end_y 
						&& k != index) {
					
					neighbours.add(population.get(k));
					
				}
				
			}	
		}
		
		return neighbours;
		
	}
	
	public static int getPoissonRandom(double mean) {
		
	    Random r = new Random();
	    double L = Math.exp(-mean);
	    int k = 0;
	    double p = 1.0;
	    do {
	        p = p * r.nextDouble();
	        k++;
	    } while (p > L);
	    return k - 1;
	}
	
	public static int[][] syngamy(int[][] g1, int[][] g2) {
		
		int sum = g1.length + g2.length;
		int[][] offspring = new int[sum][g1[0].length];
		
		for(int i = 0; i < g1.length; i++) {
			for(int j = 0; j < offspring[i].length; j++) {
				offspring[i][j] = g1[i][j];
			}
		}
		
		for(int i = g1.length; i < g2.length + g1.length; i++) {
			for(int j = 0; j < offspring[i].length; j++) {
				
				offspring[i][j] = g2[i - g1.length][j];
			}
		}
		
		return offspring;
	}
	
	public static int[] newPositions(int x, int y, int radius, int[][] map) {
		
		Random rand = new Random();
		int[] positions = new int[2];
		
		ArrayList<Integer> x_coord = new ArrayList<>();
		ArrayList<Integer> y_coord = new ArrayList<>();
		
		int start_x = x - radius, end_x = x + radius;
		int start_y = y - radius, end_y = y + radius;
		
		if(start_x < 0) {
			start_x = 199 + start_x;
		}
		
		if(end_x > 199) {
			end_x = 0 + (end_x - 199);
		}
		
		if(start_y < 0) {
			start_y = 99 + start_y;
		}
		
		if(end_y > 99) {
			end_y = 0 + (end_y - 99);
		}
		
		if(start_x < end_x && start_y < end_y) {
			
			for(int i = start_x; i <= end_x; i++) {
				x_coord.add(i);
			}
			
			for(int i = start_y; i <= end_y; i++) {
				y_coord.add(i);
			}
			
			
		}else if(start_x > end_x && start_y < end_y) {
			
			for(int i = start_x; i <= 199; i++) {
				x_coord.add(i);
			}
			
			for(int i = 0; i <= end_x; i++) {
				x_coord.add(i);
			}
			
			for(int i = start_y; i <= end_y; i++) {
				y_coord.add(i);
			}
			
		}else if(start_x < end_x && start_y > end_y) {
			
			for(int i = start_x; i <= end_x; i++) {
				x_coord.add(i);
			}
			
			for(int i = start_y; i <= 99; i++) {
				y_coord.add(i);
			}
			
			for(int i = 0; i <= end_y; i++) {
				y_coord.add(i);
			}	
			
		}else if(start_x > end_x && start_y > end_y) {
			
			for(int i = start_x; i <= 99; i++) {
				x_coord.add(i);
			}
			
			for(int i = 0; i <= end_x; i++) {
				x_coord.add(i);
			}
			
			for(int i = start_y; i <= 99; i++) {
				y_coord.add(i);
			}
			
			for(int i = 0; i <= end_y; i++) {
				y_coord.add(i);
			}
			
		}
		
		positions[0] = x_coord.get(rand.nextInt(x_coord.size()));
		positions[1] = y_coord.get(rand.nextInt(y_coord.size()));
		return positions;
	}
	
	public static double computeFitness(int[][] genotype, int[] reference) {
		
			
			double current = Double.MAX_VALUE;
			double fitness = 0.0;
		    // Iterate through each chromosome in the genome
		    for (int i = 0; i < genotype.length; i++) {
		        double count = 0.0;
		        for (int j = 0; j < genotype[i].length; j++) {
		            if (genotype[i][j] != reference[j]) {
		            	count += 1.0;
		            }
		        }
		        
		        if(count < current) {
		        	double l = Math.pow(count, 1.1);
		        	fitness = Math.exp(-0.05 * l + 1.5);
		        	
		        }
		    }
		return fitness;
	}

}
