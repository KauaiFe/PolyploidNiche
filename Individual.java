package FrogsLinear2024_3;

import java.util.Random;

public class Individual {
	
	private int x_coord;
	private int y_coord;
	private double v;
	private final double triploidFitness = 0.30;
	private final double tetraploidFitness = 0.999;
	private double mutationRate;
	private double fitness;
	private final double selectionStrength = 1.0;
	private int ploidy;
	private int matingRadius;
	private int loci;
	private int[][] genotype;
	private int[][] gametes;
	private int[] reference;
	
	public Individual(int matingRadius, int x, int y, double mutation, double v, int[][] genotype, int[] reference) {
		
		this.ploidy = genotype.length;
		this.loci = genotype[0].length;
		this.x_coord = x;
		this.y_coord = y;
		this.v = v;
		this.matingRadius = matingRadius;
		this.genotype = genotype;
		this.mutationRate = mutation;
		this.reference = reference;
		
	}
	
	private void meiosis() {
		
		this.gametes = null;
		Random rand = new Random();
		
		//Let us use a switch statement conditioned on the ploidy level
		
		switch(this.ploidy) {
		case 2:

			double u = Math.random();
			if(u < (1.0 - this.v)) {
				
				this.gametes = new int[1][this.loci]; //Haploid gamete
				
				for(int i = 0; i < this.genotype[0].length; i++) {
					
					this.gametes[0][i] = this.genotype[rand.nextInt(2)][i];
					
				}
				
			}else {
				
				this.gametes = new int[2][this.loci]; //Diploid gamete
				
				for(int i = 0; i < this.genotype[0].length; i++) {
					
					int p1 = rand.nextInt(2);
					int p2 = p1;
					
					while(p2 == p1) {
						p2 = rand.nextInt(2);
					}
					
					this.gametes[0][i] = this.genotype[p1][i];
					this.gametes[1][i] = this.genotype[p2][i];
					
				}

			}

			break;
		case 3:
			
			double phi = Math.random();

			if(phi < this.triploidFitness) {
				
				double l = Math.random();
				
				//Frequency of gamete types: haploid 25%, diploid 25% and triploid 50%
				
				if(l < 0.25) {
					
					this.gametes = new int[1][this.loci]; //Haploid gamete
					
					for(int i = 0; i < this.genotype[0].length; i++) {
						
						this.gametes[0][i] = this.genotype[rand.nextInt(3)][i];
						
					}
					
				}else if(l >= 0.25 && l < 0.5) {
					
					this.gametes = new int[2][this.loci]; //Diploid gamete
					
					for(int i = 0; i < this.genotype[0].length; i++) {
						
						int p1 = rand.nextInt(3);
						int p2 = p1;
						
						while(p2 == p1) {
							p2 = rand.nextInt(3);
						}
						
						this.gametes[0][i] = this.genotype[p1][i];
						this.gametes[1][i] = this.genotype[p2][i];
						
					}

				}else if(l >= 0.5){
					
					this.gametes = new int[3][this.loci]; //Triploid gamete
					
					for(int i = 0; i < this.genotype[0].length; i++) {
						
						this.gametes[0][i] = this.genotype[0][i];
						this.gametes[1][i] = this.genotype[1][i];
						this.gametes[2][i] = this.genotype[2][i];
						
					}
				}
			}
			
			break;
		case 4:
			
			double u2 = Math.random();
			if(u2 < this.tetraploidFitness) {
				
				int p1 = rand.nextInt(4);
				int p2 = p1;
				
				while(p2 == p1) {
					p2 = rand.nextInt(4);
				}
				
				this.gametes = new int[2][this.loci]; //Diploid gamete
				
				for(int i = 0; i < this.genotype[0].length; i++) {
					
					this.gametes[0][i] = this.genotype[p1][i];
					this.gametes[1][i] = this.genotype[p2][i];

				}
				
			}
			break;
		}
		try {
			mutation();
		}catch(NullPointerException e) {
			
		}

	}
	
	private void mutation() {
		
		for(int i = 0; i < this.gametes.length; i++) {
			for(int j = 0; j < this.gametes[0].length; j++) {
				double m = Math.random();
				if(m < this.mutationRate) {
					if(this.gametes[i][j] == 0) {
						this.gametes[i][j] = 1;
					}else {
						this.gametes[i][j] = 0;
					}
				}
			}
		}
	}
	
	public double getFitness() {
		
		return fitness;
	}
	
	public int[][] getGametes() {
		meiosis();
		return gametes;
	}

	public int getX_coord() {
		return x_coord;
	}

	public int getY_coord() {
		return y_coord;
	}

	public int getPloidy() {
		return ploidy;
	}

	public int getMatingRadius() {
		return matingRadius;
	}

	public int[][] getGenotype() {
		return genotype;
	}

}

