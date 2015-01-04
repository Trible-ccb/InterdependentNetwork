package graphs;
import java.util.HashSet;


public class NetworkNodes {

	public HashSet<Integer> neberSet;
	public int index;
	@Deprecated
	public int dependentIndex;
	public HashSet<Integer> dependentSets;
	
	public int getDegree(){
		return neberSet.size();
	}
	public NetworkNodes(){
		index = 0;
		neberSet = new HashSet<Integer>();
		dependentSets = new HashSet<Integer>();
	}
	public NetworkNodes(int i){
		this();
		index = i;
	}
	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer("" + index);
		for ( Integer i : neberSet ){
			sb.append(" " + i);
		}
		sb.append(" " + dependentSets + "\n");
		for ( Integer i : dependentSets ){
			sb.append(" " + i);
		}
		sb.append(" " + dependentIndex + "\n");
		return sb.toString();
	}
	public void copyValueBy(NetworkNodes a)
	{
		this.neberSet = new HashSet<Integer>();
		this.dependentSets = new HashSet<Integer>();
		for(Integer i : a.neberSet)
			this.neberSet.add(i);
		for(Integer i : a.dependentSets)
			this.dependentSets.add(i);
		this.dependentIndex=a.dependentIndex;
		this.index = a.index;
	}
	
}
