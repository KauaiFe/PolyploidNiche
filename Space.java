

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

/* The Artificial Genome (AG) model is implemented to initialize the program.
 * Each individual receives two copies of the AGs, thus making them diploids.
 * Each AG consists of a sequence of four nucleotides taken randomly from the set {0, 1, 2, 3}
 * Each chromosome has length chrLength;
 * A population with size 'popSize' is created.
 */

public class Space {
	
	// The length of each chromosome in the artificial genome (AG).
	private int chrLength;

	// A HashMap to store the reference genome used to build the artificial genomes for the population.
	// The key is an ID associated with specific regions in the artificial genome, and the value is the corresponding chromosome sequence.
	private HashMap<Integer, int[]> refGenome;

	// A 2D array representing the map of Australia used in the simulation.
	// Each element of the array corresponds to a specific position on the map and stores a value related to variability in average annual precipitation.
	private int[][] map;
	
	public int[] getRefGenome(int i) {
		return refGenome.get(i);
	}
	
	public int[][] getMap() {
		return map;
	}
	
	public Space(int chrLength) { // Constructor to initialise simulation
		
		this.chrLength = chrLength;
		this.map = loadMap();
		buildGenes();

	}

	/**
	 *  This method is responsible for building the artificial genomes for the map based on a reference genome.
	 *  The goal is to create genetically diverse artificial genomes for each region.
	 *  
	 *  For each cell in the map data, the code checks if the precipitation value (z) at the current cell is already present 
	 *  in the diff ArrayList. If it is found in the diff ArrayList, the boolean variable annotate is set to false, 
	 *  indicating that the current cell does not contain a unique region value. 
	 *  Otherwise, if the precipitation value (z) is not found in the diff ArrayList, annotate remains true, 
	 *  indicating that the current cell contains a unique region value.
	 *  
	 *  The mthod ensures that if the precipitation is similar between two regions on the map 
	 *  (resulting in a small difference between their unique precipitation values in the diff ArrayList), 
	 *  their corresponding chromosomes in the refGenome will be more similar. On the other hand, 
	 *  regions with significantly different precipitation values will have different chromosomes, 
	 *  leading to genetic diversity in those regions.
	 */
	private void buildGenes() {

	    // Initialize the reference genome as a HashMap with an initial capacity of 101 (estimated number of regions).
	    this.refGenome = new HashMap<>(101);

	    // Create a reference chromosome of length 'chrLength' and fill it with random nucleotides.
	    int[] chromosome = new int[this.chrLength]; // Reference chromosome
	    for (int i = 0; i < chromosome.length; i++) {
	        int nucleotide = (int) (Math.random() * 4.0);
	        chromosome[i] = nucleotide;
	    }

	    // Create an ArrayList 'diff' to store unique region values from the map data.
	    ArrayList<Integer> diff = new ArrayList<>();

	    // Loop through the map data to find unique region values and store them in the 'diff' ArrayList.
	    for (int i = 0; i < this.map.length; i++) {
	        for (int j = 0; j < this.map[i].length; j++) {
	            boolean annotate = true;
	            for (int z = 0; z < diff.size(); z++) {
	                if (this.map[i][j] == diff.get(z)) {
	                    annotate = false;
	                }
	            }
	            if (annotate) {
	                diff.add(this.map[i][j]);
	            }
	        }
	    }

	    // Sort the 'diff' ArrayList in reverse order to obtain the region values in descending order.
	    Collections.sort(diff, Collections.reverseOrder());

	    // Initialize variables to keep track of the current position in the chromosome.
	    int start = 0;
	    int[] startChrom = new int[this.chrLength];

	    // Copy the initial chromosome to 'startChrom'.
	    for (int i = 0; i < startChrom.length; i++) {
	        startChrom[i] = chromosome[i];
	    }

	    /** Loop through the 'diff' ArrayList to build the artificial genomes for each region.
	     * This  ensures that if the precipitation is similar between two compared regions on the map, 
	     * then the corresponding chromosomes for those regions will also be more similar.
	     */
	    int count = 0;
	    while (count < diff.size()) {

	        // Create a new chromosome array to store the modified chromosome for the current region.
	        int[] newChromosome = new int[this.chrLength];

	        /* The variable end is calculated as the difference between the total chromosome length (this.chrLength) 
	        * and the unique precipitation value from the diff ArrayList corresponding to the current region. */
	        int end = this.chrLength - diff.get(count);

	        // Loop through each position in the newChromosome array.
	        for (int i = 0; i < newChromosome.length; i++) {
	        	/* Within the for loop, the code checks if the current position (represented by variable i) 
            	 * is within the boundary of the current region (i.e., i >= start and i < end). 
            	 * If it is within the region boundary, it selects a random nucleotide different from the reference 
            	 * chromosome value at that position (chromosome[i]).
            	 */
	            if (i >= start && i < end) {
	            	boolean change = false;
	                while (change == false) {
	                    int nucleo = (int) (Math.random() * 4.0);
	                    if (nucleo != chromosome[i]) {
	                        newChromosome[i] = nucleo;
	                        change = true;
	                    }
	                }
	            } else {
	                // For positions outside the current region, copy the nucleotide from the reference chromosome.
	                newChromosome[i] = chromosome[i];
	            }
	        }

	        // Add the new chromosome to the reference genome with the corresponding region ID (from 'diff').
	        this.refGenome.put(diff.get(count), newChromosome);

	        // Copy the new chromosome to the original chromosome for the next iteration.
	        for (int i = 0; i < newChromosome.length; i++) {
	            chromosome[i] = newChromosome[i];
	        }

	        // Update the 'start' position for the next region.
	        start = end;
	        count++;
	    }
	}

	
	/*
	 * Load Australian map with percentages describing variability in average annual precipitation.
	 * This method reads the data from a CSV file containing precipitation percentages for each location on the map.
	 * It populates the 'map' array with the precipitation percentages for each location.
	 */
	private int[][] loadMap() {
		
		// Define the path to the CSV file containing the map data.
		String path = "/Users/feper/Documents/ThirdPaper_PHD/MapMonopolization2.csv";
		
		// Initialize variables to read data from the CSV file.
		String line = "";
		ArrayList<String> data = new ArrayList<>();
		
		// Create a 2D array to store the map data (245 rows x 176 columns).
		int[][] map = new int[200][100];
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(path));
			int count = 0;
			while((line = br.readLine()) != null) {
				data.add(line);
				String[] values = line.split(",");
				if(count != 0) {
					int x = Integer.parseInt(values[0]);
					int y = Integer.parseInt(values[1]);
					int z = Integer.parseInt(values[2]);
					map[x][y] = z;
					
				}
				count++;
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// Return the 2D 'map' array containing precipitation percentages for each location on the map.
		return map;
	}
	
}

