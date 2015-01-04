package graphs;

import java.util.HashMap;

public class Network {

	public NetworkNodes[] networkNodes;
	public int edges;
	public int nodes;
	public float avgdegree;
	public String displayname;
	public String type;
	public HashMap<Integer, Float> degreeDistribution;
	
	public Network(){
		this(0);
	}
	public Network(int size){
		nodes = size;
		init();
	}
	public void init(){
		networkNodes =  new NetworkNodes[nodes];
		for ( int i = 0 ; i < nodes ; i++ ){
			networkNodes[i] = new NetworkNodes(i);
		}
	}
	public HashMap<Integer, Float> getDegreeDistribution(){
		if (degreeDistribution == null){
			degreeDistribution = new HashMap<Integer, Float>();
			for ( int i = 0 ; i < networkNodes.length ; i++ ){
				int k = networkNodes[i].getDegree();
				Float nums = degreeDistribution.get(k);
				nums = nums == null ? 0 : nums;
				degreeDistribution.put(k, 1f / nodes + nums);
			}
//			for ( int k = 2 ; k < nodes ; k++ ){
//				float testVa = (4.0f/k/k)-(4.0f/(k+1)/(k+1));
//				degreeDistribution.put(k, testVa);
//			}
		}
		return degreeDistribution;
	}
}
