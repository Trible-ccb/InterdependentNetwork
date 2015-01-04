package graphs;



/**
 * @author Trible Chen
 *@see Newman PRL 2009
 */
public class ClusteringERGraph {

	/**
	 * @param number 
	 * @param avgdeg
	 * @param c the clustering coefficient
	 * @return
	 */
	public Network buildGraph(int number,int avgdeg,float c){
		int totaledges = avgdeg * number / 2;
		ERGraph erg = new ERGraph();
		Network er = erg.buildERGraph1(number, avgdeg);
		er.displayname = "CER_n"+number/1000+"k";
		er.type = "CER";
		//configure step
		
		
		
		System.out.println("clustered er network build over!");
		return er;
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
}
