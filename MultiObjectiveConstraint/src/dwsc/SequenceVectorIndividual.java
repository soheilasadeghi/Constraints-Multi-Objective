package wsc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import ec.EvolutionState;
import ec.Fitness;
import ec.multiobjective.MultiObjectiveFitness;

import ec.util.Parameter;
import ec.vector.VectorIndividual;

public class SequenceVectorIndividual extends VectorIndividual {

	private static final long serialVersionUID = 1L;

	private double dcost;
	private double ccost;
	private double ctime;
	private double dtime;

	private double time;
	private double cost;
	public Service[] genome;
	public List<Service> filteredgenome;

	@Override
	public Parameter defaultBase() {
		return new Parameter("sequencevectorindividual");
	}

	@Override
	/**
	 * Initializes the individual.
	 */
	public void reset(EvolutionState state, int thread) {
		WSCInitializer init = (WSCInitializer) state.initializer;
		List<Service> relevantList = init.relevantList;
		Collections.shuffle(relevantList, init.random);

		genome = new Service[relevantList.size()];
		relevantList.toArray(genome);
		this.evaluated = false;
	}

	@Override
	public boolean equals(Object ind) {
		boolean result = false;

		if (ind != null && ind instanceof SequenceVectorIndividual) {
			result = true;
			SequenceVectorIndividual other = (SequenceVectorIndividual) ind;

			for (int i = 0; i < genome.length; i++) {
				if (!genome[i].equals(other.genome[i])) {
					result = false;
					break;
				}

			}
		}
		return result;
	}

	@Override
	public int hashCode() {
		return Arrays.hashCode(genome);
	}

	@Override
	public String toString() {
		return Arrays.toString(genome);
	}

	@Override
	public SequenceVectorIndividual clone() {
    	SequenceVectorIndividual g = new SequenceVectorIndividual();
    	g.species = this.species;
    	if (this.fitness == null)
    		g.fitness = (Fitness) g.species.f_prototype.clone();
    	else
    		g.fitness = (Fitness) this.fitness.clone();
    	if (genome != null) {
    		// Shallow cloning is fine in this approach
    		g.genome = genome.clone();
    	}
    	return g;
	}

	public String toGraphString(EvolutionState state) {
		WSCInitializer init = (WSCInitializer) state.initializer;
		Graph g = createNewGraph(init.numLayers, init.startServ, init.endServ, genome, init);
		return g.toString();
	}

	public Graph createNewGraph(int numLayers, Service start, Service end, Service[] sequence, WSCInitializer init) {
		Node endNode = new Node(end);
		Node startNode = new Node(start);

        Graph graph = new Graph();
        graph.nodeMap.put(endNode.getName(), endNode);

        // Populate inputs to satisfy with end node's inputs
        List<InputNodeLayerTrio> nextInputsToSatisfy = new ArrayList<InputNodeLayerTrio>();

        for (String input : end.getInputs()){
            nextInputsToSatisfy.add( new InputNodeLayerTrio(input, end.getName(), numLayers) );
        }

        // Fulfil inputs layer by layer
        for (int currLayer = numLayers; currLayer > 0; currLayer--) {

            // Filter out the inputs from this layer that need to fulfilled
            List<InputNodeLayerTrio> inputsToSatisfy = new ArrayList<InputNodeLayerTrio>();
            for (InputNodeLayerTrio p : nextInputsToSatisfy) {
               if (p.layer == currLayer)
                   inputsToSatisfy.add( p );
            }
            nextInputsToSatisfy.removeAll( inputsToSatisfy );

            int index = 0;
            while (!inputsToSatisfy.isEmpty()){

                if (index >= sequence.length) {
                    nextInputsToSatisfy.addAll( inputsToSatisfy );
                    inputsToSatisfy.clear();
                }
                else {
                	Service nextNode = sequence[index++];
                	if (nextNode.layer < currLayer) {
	                    Node n = new Node(nextNode);
	                    //int nLayer = nextNode.layerNum;

	                    List<InputNodeLayerTrio> satisfied = getInputsSatisfiedGraphBuilding(inputsToSatisfy, n, init);

	                    if (!satisfied.isEmpty()) {
	                        if (!graph.nodeMap.containsKey( n.getName() )) {
	                            graph.nodeMap.put(n.getName(), n);
	                        }

	                        // Add edges
	                        createEdges(n, satisfied, graph);
	                        inputsToSatisfy.removeAll(satisfied);


	                        for(String input : n.getInputs()) {
	                            nextInputsToSatisfy.add( new InputNodeLayerTrio(input, n.getName(), n.getLayer()) );
	                        }
	                    }
	                }
                }
            }
        }

        // Connect start node
        graph.nodeMap.put(startNode.getName(), startNode);
        createEdges(startNode, nextInputsToSatisfy, graph);

        return graph;
    }

	public void createEdges(Node origin, List<InputNodeLayerTrio> destinations, Graph graph) {
		// Order inputs by destination
		Map<String, Set<String>> intersectMap = new HashMap<String, Set<String>>();
		for(InputNodeLayerTrio t : destinations) {
			addToIntersectMap(t.service, t.input, intersectMap);
		}

		for (Entry<String,Set<String>> entry : intersectMap.entrySet()) {
			Edge e = new Edge(entry.getValue());
			origin.getOutgoingEdgeList().add(e);
			Node destination = graph.nodeMap.get(entry.getKey());
			destination.getIncomingEdgeList().add(e);
			e.setFromNode(origin);
        	e.setToNode(destination);
        	graph.edgeList.add(e);
		}
	}

	private void addToIntersectMap(String destination, String input, Map<String, Set<String>> intersectMap) {
		Set<String> intersect = intersectMap.get(destination);
		if (intersect == null) {
			intersect = new HashSet<String>();
			intersectMap.put(destination, intersect);
		}
		intersect.add(input);
	}

	public List<InputNodeLayerTrio> getInputsSatisfiedGraphBuilding(List<InputNodeLayerTrio> inputsToSatisfy, Node n, WSCInitializer init) {
	    List<InputNodeLayerTrio> satisfied = new ArrayList<InputNodeLayerTrio>();
	    for(InputNodeLayerTrio p : inputsToSatisfy) {
            if (init.taxonomyMap.get(p.input).servicesWithOutput.contains( n.getService() ))
                satisfied.add( p );
        }
	    return satisfied;
	}

	   public void calculateSequenceFitness(int numLayers, Service end, Service[] sequence, WSCInitializer init, EvolutionState state, boolean isOperation) {
			double highestLocalTimeVariable=0.0;

	        Set<Service> solution = new HashSet<Service>();

	        filteredgenome = new ArrayList<Service>();
	        cost = 0.0;
	        ccost = 0.0;
	        dcost = 0.0;


	        List<InputTimeLayerTrio> nextInputsToSatisfy = new ArrayList<InputTimeLayerTrio>();
	        double t = end.getQos()[WSCInitializer.TIME];
	        double dt=end.getQos()[WSCInitializer.DATASETTIME];
	        for (String input : end.getInputs()){
	            nextInputsToSatisfy.add( new InputTimeLayerTrio(input, t, dt,0.0,0.0,numLayers,end.getName(),-1) );
	        }

	        // Fulfil inputs layer by layer
	        for (int currLayer = numLayers; currLayer > 0; currLayer--) {
	            // Filter out the inputs from this layer that need to fulfilled
	            List<InputTimeLayerTrio> inputsToSatisfy = new ArrayList<InputTimeLayerTrio>();
	            for (InputTimeLayerTrio p : nextInputsToSatisfy) {
	               if (p.layer == currLayer)
	                   inputsToSatisfy.add( p );
	            }
	            nextInputsToSatisfy.removeAll( inputsToSatisfy );

	            int index = 0;
	            while (!inputsToSatisfy.isEmpty()){
	                // If all nodes have been attempted, inputs must be fulfilled with start node
	                if (index >= sequence.length) {
	                    nextInputsToSatisfy.addAll(inputsToSatisfy);
	                    inputsToSatisfy.clear();
	                }
	                else {
	                Service nextNode = sequence[index++];
	                int nLayer=nextNode.layer;
	                if (nLayer < currLayer) {

	   	                List<InputTimeLayerTrio> satisfied = getInputsSatisfied(inputsToSatisfy, nextNode, init);
	   	             if (!satisfied.isEmpty()) {
                         double[] qos = nextNode.getQos();
                         if (!solution.contains( nextNode )) {
                             solution.add(nextNode);
                             filteredgenome.add(nextNode);
                             cost += qos[WSCInitializer.COST];
                             dcost += qos[WSCInitializer.DATASETCOST];
                           //  reliability *= qos[WSCInitializer.RELIABILITY];
                         }
                         highestLocalTimeVariable=0.0;
                         ccost=0.0;
                         ctime=0.0;

                         double highestLocalCostVariable=0.0;
                         double[] comarray= new double[2];

                         						for(InputTimeLayerTrio isat: satisfied)
                         {						comarray=communication.ComunicationCostandTime(isat.ID, nextNode.getID());
                         if(comarray[0]>highestLocalCostVariable) {highestLocalCostVariable=comarray[0]; ccost=highestLocalCostVariable+isat.comcost;}
                         if(comarray[1]>highestLocalTimeVariable) {highestLocalTimeVariable=comarray[1]; ctime=highestLocalTimeVariable+isat.comtime;}

                         }



                         						//now calculate distance between n and satisfieds
                         //cost between all satisfied to n
                         						//here we should consider time of links as well maximum time between satisfieds and n (max comtimes of all satisfiedstime+this max time)
                         						t = qos[WSCInitializer.TIME];
                         						dt=qos[WSCInitializer.DATASETTIME];

                         						inputsToSatisfy.removeAll(satisfied);
                         						double[] arrat=findHighestTime(satisfied);
                         						double highestT = arrat[1];
                         						double highestDT = arrat[2];
                         						for(String input : nextNode.getInputs()) {
                         							nextInputsToSatisfy.add( new InputTimeLayerTrio(input, highestT + t, highestDT+dt,ctime, ccost, nLayer, nextNode.getName(),nextNode.getID()));
                         						}
                     }
	               }
              }
          }
      }


	        // Find the highest overall time
	    			double timearray[] = findHighestTime(nextInputsToSatisfy);

	    			//particle.availability = availability;
	    			//particle.reliability = reliability;
	    			time = timearray[1];


	    			dtime=timearray[2];
	    	        if (!WSCInitializer.dynamicNormalisation || isOperation)
	    	        	finishCalculatingSequenceFitness(init, state);
	    	    }

	    	   public void finishCalculatingSequenceFitness(WSCInitializer init, EvolutionState state) {
	    		   double[] objectives = calculateFitness(cost, time, dcost, dtime, ccost,ctime, init);
	    			//init.trackFitnessPerEvaluations(f);

	    			((MultiObjectiveFitness) fitness).setObjectives(state,objectives); // XXX Move this inside the other one
	    			evaluated = true;
	    	   }

	    	   public double[] findHighestTime(List<InputTimeLayerTrio> satisfied) {
	    			double max[] = new double[3];
	    			max[0]= Double.MIN_VALUE;
	    			max[1]= 0;
	    			max[2]= 0;
	    			//here we should modify and add comtime
	    			for (InputTimeLayerTrio p : satisfied) {
	    				if ((p.time+p.datatime) > max[0])
	    					{max[0] = p.time+p.datatime;
	    					max[1]=p.time;
	    					max[2]=p.datatime;

	    					}
	    			}

	    			return max;
	    		}


		public double[] calculateFitness(double c, double t, double dc, double dt,double cc, double ct, WSCInitializer init) {

	        double[] objectives = new double[2];
	        //objectives[GraphInitializer.AVAILABILITY] = a;
	        //objectives[GraphInitializer.RELIABILITY] = r;
	        //objectives[WSCInitializer.TIME] = t;
	        //objectives[WSCInitializer.COST] = c;
	        dt = normaliseDTime(dt,init);
	        dc = normaliseDCost(dc,init);
//System.out.println(dt);
	        ct = normaliseCTime(ct,init);
	        cc = normaliseCCost(cc,init);


        t = normaliseTime(t, init);
        c = normaliseCost(c, init);
       // objectives[0]=(dt+ct+t)/3.0;
       // objectives[1]=(dc+cc+c)/3.0;
        objectives[0]=t;
         objectives[1]=c;
	        return objectives;
		}

		private double normaliseCCost(double cc, WSCInitializer init) {
			if (init.maxCCost - init.minCCost == 0.0)
				return 1.0;
			else
				return (init.maxCCost-cc)/(init.maxCCost  - init.minCCost);
		}

		private double normaliseCTime(double ct, WSCInitializer init) {
			if (init.maxCTime - init.minCTime == 0.0)
				return 1.0;
			else
				return (init.maxCTime-ct)/(init.maxCTime  - init.minCTime);
		}

		private double normaliseDCost(double dc, WSCInitializer init) {
			if (init.maxDCost - init.minDCost == 0.0)
				return 1.0;
			else
				return (init.maxDCost - dc)/(init.maxDCost - init.minDCost);
		}

		private double normaliseDTime(double dt, WSCInitializer init) {
			if (init.maxDTime - init.minDTime == 0.0)
				return 1.0;
			else
				return (init.maxDTime - dt)/(init.maxDTime - init.minDTime);
		}

	/*	private double normaliseAvailability(double availability, WSCInitializer init) {
			if (init.maxAvailability - init.minAvailability == 0.0)
				return 1.0;
			else
				return (availability - init.minAvailability)/(init.maxAvailability - init.minAvailability);
		}

		private double normaliseReliability(double reliability, WSCInitializer init) {
			if (init.maxReliability- init.minReliability == 0.0)
				return 1.0;
			else
				return (reliability - init.minReliability)/(init.maxReliability - init.minReliability);
		}
*/
		private double normaliseTime(double time, WSCInitializer init) {
			if (init.maxTime - init.minTime == 0.0)
				return 1.0;
			else
				return (init.maxTime - time)/(init.maxTime - init.minTime);
		}

		private double normaliseCost(double cost, WSCInitializer init) {
			if (init.maxCost - init.minCost == 0.0)
				return 1.0;
			else
				return (init.maxCost - cost)/(init.maxCost - init.minCost);
		}

		public List<InputTimeLayerTrio> getInputsSatisfied(List<InputTimeLayerTrio> inputsToSatisfy, Service n, WSCInitializer init) {
		    List<InputTimeLayerTrio> satisfied = new ArrayList<InputTimeLayerTrio>();
		    for(InputTimeLayerTrio p : inputsToSatisfy) {
	            if (init.taxonomyMap.get(p.input).servicesWithOutput.contains( n ))
	                satisfied.add( p );
	        }
		    return satisfied;
		}

		public void setDTime(double dtime) {
			this.dtime = dtime;
		}

		public void setDCost(double dcost) {
			this.dcost = dcost;
		}
		public void setCTime(double ctime) {
			this.ctime = ctime;
		}

		public void setCCost(double ccost) {
			this.ccost = ccost;
		}
		public void setTime(double time) {
			this.time = time;
		}

		public void setCost(double cost) {
			this.cost = cost;
		}

		public double getCTime() {
			return ctime;
		}

		public double getCCost() {
			return ccost;
		}
		public double getDTime() {
			return dtime;
		}

		public double getDCost() {
			return dcost;
		}

		public double getTime() {
			return time;
		}

		public double getCost() {
			return cost;
		}

}
