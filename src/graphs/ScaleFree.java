package graphs;
import java.util.Date;
import java.util.Random;

/**
 * @author Trible Chen
 * see price model
 */
@Deprecated
public class ScaleFree {
	
	/**
	 * @param nodeNumber the total number of nodes
	 * @param m the number of edges adding every time
	 * @return a network by nodes
	 */
	public NetworkNodes[] BuildSFNetwork(int nodeNumber,int m)
	{
		NetworkNodes[] Networks = new NetworkNodes[nodeNumber];
		
		NetworkNodes tempNode=new  NetworkNodes(0);
		Networks[0]=tempNode;
		tempNode=new  NetworkNodes(1);
		Networks[1]=tempNode;
		
		tempNode=new  NetworkNodes(2);
		Networks[2]=tempNode;
		
		AddOneLine(0,1,Networks);
		AddOneLine(1,0,Networks);
		
		AddOneLine(0,2,Networks);
		AddOneLine(2,0,Networks);
		
		AddOneLine(1,2,Networks);
		AddOneLine(2,1,Networks);
		
		for(int i=3;i<nodeNumber;i++)
		{
			
			tempNode=new  NetworkNodes(i);
			Networks[i]=tempNode;
			
			int index1=GetANodeByProbabliy(i,m,Networks);
			int index2=GetANodeByProbabliy(i,m,Networks);
			while(index1==index2)
			{
				index1=GetANodeByProbabliy(i,m,Networks);
				index2=GetANodeByProbabliy(i,m,Networks);
			}
			
			AddOneLine(i,index1,Networks);
			AddOneLine(index1,i,Networks);
			
			AddOneLine(i,index2,Networks);
			AddOneLine(index2,i,Networks);
			
		}
//		for(int i=0;i<nodeNumber;i++)
//		{
//			Networks[i].SetItParent(Networks);
//		}
		System.out.println("sf network build over!");
		return Networks;
	}
	void AddOneLine(int i,int j , NetworkNodes[] net){
		net[i].neberSet.add(j);
	}
	private int GetANodeByProbabliy(int i,int m,NetworkNodes[] Networks)
	{
//		Random rand=new Random(new Date().getTime());
		double pro = Math.random();
		float totald = (float)((3+(i-3)*2)*2);
		float nodesize=(float)(Networks[0].neberSet.size())/totald;
		int j=0;
		while(pro>nodesize)
		{
			j++;
			if(j>=i-1)break;
			nodesize+=(float)Networks[j].neberSet.size()/totald;
		}
		
		return j;
	}

//	public int FindItOldParent(NetworkNodes[] Networks,NetworkNodes node)
//	{
//		int OldParentId=node.parentId;
//		while(OldParentId!=node.index)
//		{
//			node=Networks[OldParentId];
//			OldParentId=node.parentId;
//		}
//		return OldParentId;
//	}
}
