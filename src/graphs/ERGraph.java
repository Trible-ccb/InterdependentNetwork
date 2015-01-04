package graphs;
import java.util.HashMap;
import java.util.HashSet;



public class ERGraph extends Network{

	/**
	 * @param n the network's size
	 * @param avgdeg the average degree of the network
	 * @param p the probability of connecting two nodes
	 * @return network
	 */
	public ERGraph(){
		super();
	}
	public ERGraph(int s){
		super(s);
	}
	@Deprecated
	public Network buildERgraph2(int n,int avgdeg){
		float p = (float)avgdeg / (n - 1);
		if ( n < 2 || avgdeg < 0 || p < 0 )throw new IllegalArgumentException();
		Network er = new ERGraph(n);
		NetworkNodes[] Networks = er.networkNodes;
		int cntedge = 0;
		for ( int i = 0 ; i < n ; i++ ){
			for ( int j = 0 ; j < n / 2 ; j++ ){
				if ( i !=j ){
					double pro = Math.random();
					if ( pro < p ){
						cntedge++;
						AddOneLine(i, j, Networks);
						AddOneLine(j, i, Networks);
					}
//					if ( cntedge * 2 == avgdeg * n ){
//						er.edges = cntedge;
//						er.avgdegree = avgdeg;
//						return er;
//					}
				}
			}
		}
		er.edges = cntedge;
		er.avgdegree = cntedge * 2 / (float)n;
		return er;
	}
	/**
	 * @param number 
	 * @param avgdeg
	 * @return
	 */
	public Network buildERGraph1(int number,int avgdeg){
		this.nodes = number;
		init();
		int totaledges = avgdeg * number / 2;
		NetworkNodes[] Networks = networkNodes;
		displayname = "E_n"+number/1000+"k" + "ad"+avgdeg;
		type = "ER";
		for(int i=0;i<number;i++)
		{
			Networks[i]=new NetworkNodes(i);
		}
		for(int i=0;i<totaledges;i++)
		{
			int index1=0;
			int index2=index1;
			while(index1==index2||Networks[index1].neberSet.contains(index2))//随机选两个不同的点，连接两条边
			{
				index1=GetANodeByProbabliy(number);
				index2=GetANodeByProbabliy(number);
				
			}
			AddOneLine(index1, index2, Networks);
			AddOneLine(index2, index1, Networks);
		}
		avgdegree = avgdeg;
		System.out.println("er network build over!");
		return this;
	}
	
	private void AddOneLine(int i,int j , NetworkNodes[] net){
		net[i].neberSet.add(j);
	}
	private void removeOneLine(int i,int j , NetworkNodes[] net){
		net[i].neberSet.remove(j);
	}
	public static int GetANodeByProbabliy( int n )
	{
		int edges = n;
		double index = Math.random() * edges;
		return ((int) index % edges);
	}
	@Override
	public HashMap<Integer, Float> getDegreeDistribution() {
//		if (degreeDistribution == null){
//			degreeDistribution = new HashMap<Integer, Float>();
//			for ( int k = 0 ; k < nodes ; k++ ){
//				float avg = this.avgdegree;
//				double mutik = 1;
//				double testVa;
//				if ( k > avg*25 ){
//					testVa = 0;
//				} else {
//					for ( int i = 1 ; i <= k ; i++ ){
//						mutik *= i;
//					}
//					testVa = Math.pow(avg, k) / mutik * Math.pow(Math.E, -avg);
//					degreeDistribution.put(k, (float)testVa);
//				}
//			}
//		}
//		return degreeDistribution;
		return super.getDegreeDistribution();
	}
}
