package graphs;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Trible Chen
 *无标度网络的一般模型
 */
public class Price extends Network{

	float lamda = 0;
	/**
	 * @param n the size of the network
	 * @param m here i define m0 = m
	 * @param a
	 * 			p = m / (m + a);lamada = 2 + a/m
	 * @return
	 */
	public Price(){
		super();
	}
	
	public Price(int s){
		super(s);
	}
//	public int size,m;
//	public float lamada;
//	public NetworkNodes[] networkNodes;
	
	public Network buildNetwork(int n,int m,float a){
		this.nodes = n;
		init();
		NetworkNodes[] networks = networkNodes;
		int m0 = m;
		float p = m / ( a + m );
		this.lamda = 2 + a / m;
		this.displayname = "sf"+lamda+"_n"+n/1000+"k"+"ad"+2*m0;
		this.avgdegree = m * 2;
		this.type = "SF";
//		size = n;
//		this.m = m;
//		this.lamada = lamda;
//		the list contains the node index,if the node has two degree ,the index of the node will be add twice.
		List<Integer> randomIndex = new ArrayList<Integer>();
//		initial with a m0 circle network
		for ( int i = 0 ; i < m0 ; i++ ){
			networks[i] = new NetworkNodes(i);
		}
		for ( int i = 0 ; i < m0 ; i++ ){
			int nebindx = (i + 1 + m0) % m0;
			networks[i].neberSet.add((i+1)%m0);
			networks[nebindx].neberSet.add(i);
			randomIndex.add(i);
			randomIndex.add(nebindx);
		}
		for ( int i = m0 ; i < n ; i++ ){//for each new node
			NetworkNodes newnode = new NetworkNodes(i);
			networks[i] = newnode;
			Set<Integer> hadselect = new HashSet<Integer>();
			for ( int j = 0 ; j < m ; ){//choose m nodes for the new node
				double r = Math.random();
				int n1;
				if ( r < p ){//choose node from the random index list
					n1 = (int) Math.floor(Math.random() * (randomIndex.size() - 1));
					n1 = randomIndex.get(n1);
				} else {//choose node from the network list
					n1 = (int) Math.floor(Math.random() * (i - 1));
				}
				if ( hadselect.contains(n1) ){
					continue;
				} else {
					j++;
					hadselect.add(n1);
				}
			}
			newnode.neberSet.addAll(hadselect);
			for ( Integer nebindx : hadselect ){
				networks[nebindx].neberSet.add(i);
				randomIndex.add(nebindx);
			}

		}
		return this;
	}
	@Override
	public HashMap<Integer, Float> getDegreeDistribution() {
//		if (degreeDistribution == null){
//			degreeDistribution = new HashMap<Integer, Float>();
//			for ( int k = 2 ; k < nodes ; k++ ){
//				float testVa = ((lamda - 1)*(lamda - 1)/k/k)
//						-((lamda - 1)*(lamda - 1)/(k+1)/(k+1));
//				degreeDistribution.put(k, testVa);
//			}
//		}
//		return degreeDistribution;
		return super.getDegreeDistribution();
	}
}
