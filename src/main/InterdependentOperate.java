package main;
import graphs.ERGraph;
import graphs.Network;
import graphs.NetworkNodes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.Set;
import java.util.concurrent.LinkedBlockingQueue;

import javax.print.attribute.standard.Finishings;




public class InterdependentOperate {
	static double eps = 0.001;
	
	private static boolean isEqual(double p,double p1){
		if ( Math.abs(p-p1)<eps){
			return true;
		}
		return false;
	}
	public static void dependentTwoNetRandom(NetworkNodes[] A,NetworkNodes[] B){
		int NetworkSize = A.length;
		if ( NetworkSize != B.length ){
			return;
		}
		ArrayList<Integer> indexs = new ArrayList<Integer>();
		for ( int i = 0 ; i < NetworkSize ; i++ ){
			NetworkNodes node = B[i];
			indexs.add(node.index);
		}
//		boolean[] isused = new boolean[NetworkSize];
		for(int i=0;i<NetworkSize;i++)
		{
			int j=ERGraph.GetANodeByProbabliy(indexs.size());
//			while(isused[j])
//			{
//				j=ERGraph.GetANodeByProbabliy(NetworkSize);
//			}
//			j = B[indexs.get(j)].index;
//			isused[j]=true;
			j = indexs.remove(j);
			A[i].dependentIndex=j;
			B[j].dependentIndex=i;
		}
		System.out.println("Init dependent over!");
	}
//	φ(k)=k^α*∑kj^β
	static double fix = 1;
	public static HashMap<Integer, Double> funFAI(Network net,Network net2,float alpha,float beta){
		
		HashMap<Integer, Double> result = new HashMap<Integer, Double>();
		HashMap<Integer, Integer> cnt = new HashMap<Integer, Integer>();
		HashMap<Integer, HashMap<Integer, Double>> jointP = new HashMap<Integer, HashMap<Integer, Double>>();
		
		int size = net.nodes;
		double sum = 0;
		for (int i = 0 ;i < size ; i++){
			int k1 = net.networkNodes[i].getDegree();
			int k2 = net2.networkNodes[net.networkNodes[i].dependentIndex].getDegree();
			HashMap<Integer, Double> k1map = jointP.get(k1);
			if ( k1map == null )k1map = new HashMap<Integer, Double>();
			jointP.put(k1, k1map);
			Double k12num = k1map.get(k2);
			k12num = k12num == null ? 1 : k12num+1;
			k1map.put(k2, k12num);
//			Double va = result.get(k1);
//			double ret = Math.pow(k1+1, alpha)*Math.pow(k2+1, beta);
//			sum += ret;
//			va = va == null ? 0 : va;
//			result.put(k1, va+ret);
//			Integer ia = cnt.get(k1);
//			ia = ia == null ? 0 : ia;
//			cnt.put(k1, 1+ia);
		}
		for ( Integer k1 : jointP.keySet() ){
			HashMap<Integer, Double> k1v = jointP.get(k1);
			double k1sum = 0;
			for ( Entry<Integer, Double> k12num:k1v.entrySet() ){
				Integer k2 = k12num.getKey();
				Double numP = k12num.getValue() / size * Math.pow(k2+fix, beta);
				k1sum+=numP;
			}
			result.put(k1, k1sum*Math.pow(k1+fix, alpha));
		}
//		HashMap<Integer, Float> distr = net.getDegreeDistribution();
//		for ( Entry<Integer, Float> entry : distr.entrySet()) {
//			int k = entry.getKey();
//			if( k == 0 ){
//				result.put(k,0.0+fix);
//			} else {
//				result.put(k,Math.pow(k+fix, alpha));
//			}
//		}
		return result;
	}
//	变量t
	private static double getVariableT(HashMap<Integer, Double> faiKs,Network net,Network net2,float p){
		HashMap<Integer, Float> distr = net.getDegreeDistribution();
		double x = 0.5;
		double start = 0.00001;
		double end = 1;
		
		while(true){
			//二分方法查找t
			x = (start + end) / 2;
			if ( isEqual(x, 0) ){
				break;
			}
			double sum = 0;
//			用于判断单调性
			double sumStart = 0;
			double sumEnd = 0;
			for ( Entry<Integer, Float> entry : distr.entrySet() ){
				int k = entry.getKey();
				double fai = faiKs.get(k);
				sum += entry.getValue()*Math.pow(x, fai);
				sumStart += entry.getValue()*Math.pow(start, fai);
				sumEnd += entry.getValue()*Math.pow(end, fai);
			}
			
			if ( Math.abs(p-sum)<0.05 ){
				break;
			} else {
				if( sumStart < sumEnd ){//单调递增
					if( sum > p ){
						end = x;
					} else {
						start = x;
					}
				} else {//单调递减
					if( sum > p ){
						start = x;
					} else {
						end = x;
					}
				}

			}
		}
		double t1 = Math.log(p)/net.avgdegree + 1;
		return fixRet(x);
	}
//	Gib(x)
	private static double funGb(double x,double t,HashMap<Integer, Double> faiKs,Network net,Network net2,float p){
		HashMap<Integer, Float> distr = net.getDegreeDistribution();
		double sum = 0;
		for ( Entry<Integer, Float> entry : distr.entrySet() ){
			int k = entry.getKey();
			double fai = faiKs.get(k);
			sum += entry.getValue()
					*Math.pow(t, fai)
					*Math.pow(x, k);
		}
		if ((sum /p)>1){
			int a = 0;
		}
		return fixRet(sum /p);
	}
//	Gi1(x)
	@SuppressWarnings("unused")
	private static double funGone(double x,double t,double pyp,HashMap<Integer, Double> faiK,Network net,Network net2,float p){
		double fenzi,fenmu;
//		double pyp = getProbabilityPyipiao(t,net,net2,alpha,p,beta);
		double nx = 1 + pyp / p * ( x - 1 );
		fenzi = derivativeFunGb(nx,t, faiK,net,net2, p);
		fenmu = derivativeFunGb(1,t, faiK,net,net2, p);
		if ( fenzi / fenmu>1){
			int a =0 ;
		}
		return fenzi / fenmu;
	}
//	Gib(x) 导函数
	private static double derivativeFunGb(double x,double t,HashMap<Integer, Double> faiKs,Network net,Network net2,float p){
		HashMap<Integer, Float> distr = net.getDegreeDistribution();
		double sum = 0;
		for ( Entry<Integer, Float> entry : distr.entrySet() ){
			int k = entry.getKey();
			double fai = faiKs.get(k);
			sum += k*entry.getValue()
					*Math.pow(t, fai)
					*Math.pow(x,k-1);
		}
		
		return sum / p;
	}
//	~p 变量p一飘
	private static double getProbabilityPyipiao(double t,HashMap<Integer, Double> faiKs,Network net,Network net2,float p){
		HashMap<Integer, Float> distr = net.getDegreeDistribution();
		double sumFenmu = net.avgdegree;
		double sumFenzi = 0;
		for ( Entry<Integer, Float> entry : distr.entrySet() ){
			int k = entry.getKey();
			double fai = faiKs.get(k);
			sumFenzi += k*entry.getValue()
					*Math.pow(t, fai);
		}
		if ( sumFenzi / sumFenmu>1){
			int a =0 ;
		}
		return sumFenzi / sumFenmu;
	}
	private static double fixRet(double ret){
		if ( ret>1 ){
			ret = 1;
		} else if ( ret < 0 ){
			ret = 0;
		}
		return ret;
	}
//	gi(p)
	private static double funGp(double x,double t,HashMap<Integer, Double> faiKs,Network net,Network net2,float p){
		double startFA = x;
		double newFA = 0;
		double pyp = getProbabilityPyipiao(t,faiKs,net,net2,p);
		while ( true ){
			newFA = funGone(1-x*(1-startFA),t,pyp,faiKs,net,net2,p);
			if ( isEqual(newFA, startFA) ){
				break;
			} else {
				startFA = newFA;
			}
		}
		double nx = pyp/p*x*(newFA - 1)+1;
		double ret = 1-funGb(nx,t,faiKs, net,net2, p);
		return ret;
	}
//	临界值Xc
	private static double getXc(HashMap<Integer, Double> faiKs,HashMap<Integer, Double> faiKs2,Network net,Network net2,float p){
		double startX = p;
		double newX = 0;
		double t1 = getVariableT(faiKs, net, net2, p);
		double t2 = getVariableT(faiKs2, net2, net, p);
		while ( true ){
			double y = p * funGp(startX,t2,faiKs2,net2,net,p);
			y = fixRet(y);
			newX = p * funGp(y,t1,faiKs,net,net2,p);
			newX = fixRet(newX);
			if (isEqual(newX, startX)){
				break;
			} else {
				startX = newX;
			}
		}
		return fixRet(newX);
	}
	@Deprecated
	private static double getPc(Network net,Network net2,
			float alpha,float p,float beta,double xc){
		double y = 0;
		double np = p;
		while ( true ){
//			y = np*funGp(xc, net, net2, alpha, p, beta);
			break;
		}
		return np;
	}
//	p∞
	private static double getPinfinite(HashMap<Integer, Double> faiKs2,Network net,Network net2,float p,double xc){
//		double pi = 0.5;
//		double t = getVariableT(net, net2, alpha, p, beta);
//		while ( true ){
//			double tmp = (1-Math.pow(Math.E, -net.avgdegree*t*t*pi));
//			double newpi = p*p*(tmp*tmp);
//			if ( isEqual(newpi, pi) ){
//				break;
//			} else {
//				pi = newpi;
//			}
//		}
//		return pi;
		double t2 = getVariableT(faiKs2, net2, net, p);
		return fixRet(xc*funGp(xc,t2,faiKs2, net2, net, p));
	}
//	级联失效理论结果
	public static double cascadeDependentNetInTheory(
			HashMap<Integer, Double> faiKs,HashMap<Integer, Double> faiKs2,Network A,Network B,float initialRemove){
		float p = 1 - initialRemove;
		double xc = getXc(faiKs,faiKs2, A,B, p);
		double pinfinite = getPinfinite(faiKs2, A,B, p, xc);
		return pinfinite;
	}
	
	
	
	
	
//	级联恢复过程
	/**
	 * @param A
	 * @param B
	 * @param initialRecover 初始恢复顶点比例
	 * @param recoverPro 链接边恢复概率
	 */
	public static void cascadeRecovery(
			NetworkNodes[] A,NetworkNodes[] B,float initialRecover,float recoverPro){
		List<NetworkNodes> networkNodes = new ArrayList<NetworkNodes>(Arrays.asList(A));
		Set<Integer> recoeringIndex = new HashSet<Integer>();
		Queue<Integer> recoveryNodesIndexQueue = new LinkedBlockingQueue<Integer>();
		int len = A.length;
		int iniRecoveryNum = (int) (initialRecover * len);
		for ( int i = 0 ; i < iniRecoveryNum ; i++){//初始化恢复一定比例的A网络
			int remainSize = networkNodes.size();
			int chooseidx = ((int) (Math.random() * remainSize) % remainSize);
			NetworkNodes node = networkNodes.remove(chooseidx);
			recoveryNodesIndexQueue.add(node.index+1);
			recoeringIndex.add(node.index+1);
		}
		Set<Integer> visitedNodes = new HashSet<Integer>();
		while (!recoveryNodesIndexQueue.isEmpty()){
			int index = recoveryNodesIndexQueue.poll();
			if ( visitedNodes.contains(index) )continue;
			visitedNodes.add(index);
			int fixindex = 1;
			NetworkNodes[] nodes,nodes2;
			if ( index > 0 ){//is A
				nodes = A;
				nodes2 = B;
				index = index - 1;
				fixindex = -1;
			} else {
				fixindex = 1;
				nodes = B;
				nodes2 = A;
				index = -(index+1);
			}
			NetworkNodes node = nodes[index];
//			恢复依赖节点
			NetworkNodes depandeNode = nodes2[node.dependentIndex];
			recoveryNodesIndexQueue.add(depandeNode.index*fixindex+fixindex);
			recoeringIndex.add(depandeNode.index*fixindex+fixindex);
			Set<Integer> neb = node.neberSet;
			Set<Integer> needremoved = new HashSet<Integer>();
			for ( Integer nebIdx : neb){
				if ( !recoeringIndex.contains(-nebIdx*fixindex-fixindex) ){
//					按概率恢复邻居节点
					double randomrecoverypro = Math.random();
					boolean canrecovery = randomrecoverypro < recoverPro;
					if ( !canrecovery ){
						needremoved.add(nebIdx);
					} else {
						recoveryNodesIndexQueue.add(-nebIdx*fixindex-fixindex);
						recoeringIndex.add(-nebIdx*fixindex-fixindex);
					}
				}
			}
			for (Integer nebidx2 : needremoved ){
				nodes[nebidx2].neberSet.remove(node.index);
			}
			neb.removeAll(needremoved);
		}
		for ( int i = 0 ; i < len ; i++ ){
			int aidx = i+1;
			if ( !visitedNodes.contains(aidx)){
				NetworkNodes nodea = A[i];
				for ( Integer idx : nodea.neberSet ){
					A[idx].neberSet.remove(i);
				}
				nodea.neberSet.clear();
				nodea.dependentIndex = -1;
			}
			int bidx = -(i+1);
			if ( !visitedNodes.contains(bidx)){ 
				NetworkNodes nodeb = B[i];
				for ( Integer idx : nodeb.neberSet ){
					B[idx].neberSet.remove(i);
				}
				nodeb.neberSet.clear();
				nodeb.dependentIndex = -1;
			}
		}
		HashMap<Integer, NetworkNodes> Amap = new HashMap<Integer, NetworkNodes>();
		HashMap<Integer, NetworkNodes> Bmap = new HashMap<Integer, NetworkNodes>();
		for ( int i = 0 ; i < len ; i++ ){
			if ( A[i].dependentIndex != -1 ){
				Amap.put(i, A[i]);
			}
			if ( B[i].dependentIndex != -1 ){
				Bmap.put(i, B[i]);
			}
		}
		boolean canremove = true;
		while ( canremove ){
			int r1 = removeEdgeNotInTheSameCluster(Amap, Bmap);
			int r2 = removeEdgeNotInTheSameCluster(Bmap, Amap);
			canremove = !(0 == r2 && r1 ==0); 
		}
		canremove = false;
	}
	
	
	
	
//	级联失效模拟过程
	public static void cascadeDependentNet(
			NetworkNodes[] A,NetworkNodes[] B,float alpha,float initialRemove,float beta){
		
		Set<NetworkNodes> tmpA = new HashSet<NetworkNodes>();
		int initialremoveNum = (int) (initialRemove*A.length);
		//正数表示A中的节点 ，负数表示B中的点
		Queue<Integer> initialremoveIndex = new LinkedBlockingQueue<Integer>(initialremoveNum*2);
		
		double fenmu=0;
		for(int j=0;j<A.length;j++)
		{
			tmpA.add(A[j]);
			
			if ( fix == 0 && (A[j].neberSet.size()==0
					|| B[A[j].dependentIndex].neberSet.size()==0)){
				continue;
			}
			double jPro = Math.pow(A[j].neberSet.size()+fix, alpha)
						*Math.pow(B[A[j].dependentIndex].neberSet.size()+fix,beta);
			fenmu+=jPro;
		}

		for(int i=0;i<initialremoveNum;)
		{
			double pros=Math.random();
			double thisPro = 0;
			for(NetworkNodes e : tmpA)
			{
				int j = e.index;
				double jPro = 0;
				if ( fix == 0 && (A[j].neberSet.size()==0
						|| B[A[j].dependentIndex].neberSet.size()==0)){
					continue;
				}
				jPro = (Math.pow(A[j].neberSet.size()+fix, alpha)
						*Math.pow(B[A[j].dependentIndex].neberSet.size()+fix,beta));
				thisPro+=jPro / fenmu;
				if(thisPro>pros)
				{
					i++;
					tmpA.remove(e);
					fenmu -= jPro;
					initialremoveIndex.add(A[j].index + 1);//avoid the zero in both A and B
					initialremoveIndex.add(-B[A[j].dependentIndex].index - 1);
					break;
				}
			}
		}
		while ( !initialremoveIndex.isEmpty() ){
			int index = initialremoveIndex.poll();
			NetworkNodes node[],node2[];
			boolean isANetwork;
			if ( index > 0 ){//retrieve the real index
				index--;
				node = A;
				node2 = B;
				isANetwork = true;
			} else {
				index = (-index) - 1 ;
				isANetwork = false;
				node2 = A;
				node = B;
			}
			for ( Integer idx : node[index].neberSet ){
				node[idx].neberSet.remove(index);
			}
			node[index].neberSet.clear();//remove connect link
			

//			if ( node[index].dependentIndex != -1 ){//remove dependent link
//				if ( isANetwork ){//then the B network node need to add
//					initialremoveIndex.add( -(node[index].dependentIndex + 1) );
//				} else {
//					initialremoveIndex.add(node[index].dependentIndex + 1);
//				}
//				node2[node[index].dependentIndex].dependentIndex = -1;
//			}
			
			node[index].dependentIndex = -1;
		}
//		HashMap<Integer, NetworkNodes> Amap = new HashMap<Integer, NetworkNodes>();
//		HashMap<Integer, NetworkNodes> Bmap = new HashMap<Integer, NetworkNodes>();
//		for ( int i = 0 ; i < A.length ; i++ ){
//			if ( A[i].dependentIndex != -1 ){
//				Amap.put(i, A[i]);
//			}
//			if ( B[i].dependentIndex != -1 ){
//				Bmap.put(i, B[i]);
//			}
//		}
		boolean canremove = true;
		while ( canremove ){
			int r1 = removeEdgeWhenNotInTheSameCluster(A, B);
			int r2 = removeEdgeWhenNotInTheSameCluster(B, A);
//			int r1 = removeEdgeNotInTheSameCluster(Amap, Bmap);
//			int r2 = removeEdgeNotInTheSameCluster(Bmap, Amap);
			canremove = !(0 == r2 && r1 ==0); 
		}

	}
	private static int removeEdgeNotInTheSameCluster(HashMap<Integer, NetworkNodes> A,HashMap<Integer, NetworkNodes> B){
		//remove links of B which connect different clusters in A
		int removecount = 0;
		Set<Integer> cluster = null;
		for ( Iterator<Entry<Integer, NetworkNodes>> iter = B.entrySet().iterator();iter.hasNext();){
			Entry<Integer, NetworkNodes> entry = iter.next();
			NetworkNodes b = entry.getValue();
			int dependidx1 = b.dependentIndex;
			if (dependidx1 == -1){
				iter.remove();
				continue;
			}
			Set<Integer> needremove = new HashSet<Integer>();
			if ( cluster == null || !cluster.contains(dependidx1)){
				cluster = getComponentByIndex(A,dependidx1);
			}
			for ( Integer bidx : b.neberSet ){
				if ( B.get(bidx) == null 
						|| !cluster.contains(B.get(bidx).dependentIndex)){
					needremove.add(bidx);
					if ( B.get(bidx) != null  )
					B.get(bidx).neberSet.remove(b.index);
				}
			}
			removecount += needremove.size();
			b.neberSet.removeAll(needremove);
			if ( b.neberSet.isEmpty() ){
				int aidx = b.dependentIndex;
				A.get(aidx).neberSet.clear();
				A.get(aidx).dependentIndex = -1;
				b.dependentIndex = -1;
				A.remove(aidx);
				iter.remove();
			}
		}
		return removecount;
	}
	private static int removeEdgeWhenNotInTheSameCluster(NetworkNodes[] A,NetworkNodes[] B){
		//remove links which connect different clusters
		int removecount = 0;
		boolean[] isTestB = new boolean[B.length];
		
		for ( int i = 0 ; i < B.length; i++ ){
			if(isTestB[i])continue;
			
			isTestB[i]=true;
			
			NetworkNodes b = B[i];
			int dependidx1 = b.dependentIndex;
			
			if ( dependidx1 == -1)continue;//already remove
			
			//NetworkNodes a1 = A[dependidx1];
			
			Boolean[] isConnected = BFS(dependidx1,A);
			
			for(int j=0;j<isConnected.length;j++)
			{
				if(isConnected[j]==false)continue;
				if(A[j].dependentIndex==-1)continue;
				
				isTestB[A[j].dependentIndex]=true;
				b = B[A[j].dependentIndex];
				
				
				Set<Integer> needremove = new HashSet<Integer>();
				for ( Integer bidx : b.neberSet ){
					
					int dependidx2 = B[bidx].dependentIndex;
					if ( dependidx2 != -1 && isTestB[bidx] ==false ){
						NetworkNodes a2 = A[dependidx2];
						if(isConnected[a2.index])
						{
						}
						else {
							B[bidx].neberSet.remove(b.index);
							needremove.add(bidx);
						}
					}
					
				}
				removecount += needremove.size();
				b.neberSet.removeAll(needremove);
			}
		}
		return removecount;
	}
	private static Boolean[] BFS(int index,NetworkNodes[] nodes)
	{
		Queue<Integer> findneb = new LinkedBlockingQueue<Integer>();
		
		Boolean[] isused = new Boolean[nodes.length];
		for(int i=0;i<isused.length;i++)
			isused[i]=new Boolean(false);
		findneb.add(index);
		isused[index] = true;
		while( !findneb.isEmpty() ){
			int nebidx = findneb.poll();
			for ( Integer j : nodes[nebidx].neberSet ){
				if (!isused[j]){
					findneb.add(j);
					isused[j] = true;
				}
			}
		}
		return isused;
	}
	public static Set<Integer> getComponentByIndex(HashMap<Integer, NetworkNodes> nodes,int start){
		Set<Integer> set = new HashSet<Integer>();
		Queue<Integer> findneb = new LinkedBlockingQueue<Integer>();
		findneb.add(start);
		set.add(start);
		while( !findneb.isEmpty() ){
			int nebidx = findneb.poll();
			if ( nodes.get(nebidx) == null || nodes.get(nebidx).neberSet == null )continue;
			for ( Integer j : nodes.get(nebidx).neberSet ){
				if (!set.contains(j.intValue())){
					findneb.add(j);
					set.add(j);
				}
			}
		}
		return set;
	}
	public static Set<Integer> getComponentByIndex(NetworkNodes[] nodes,int start){
		Set<Integer> set = new HashSet<Integer>();
		Queue<Integer> findneb = new LinkedBlockingQueue<Integer>();
		findneb.add(start);
		set.add(start);
		while( !findneb.isEmpty() ){
			int nebidx = findneb.poll();
			for ( Integer j : nodes[nebidx].neberSet ){
				if (!set.contains(j.intValue())){
					findneb.add(j);
					set.add(j);
				}
			}
		}
		return set;
	}
	
	private static boolean isConnectedInNetwork(int index1,int index2,NetworkNodes[] nodes){
		NetworkNodes n1 = nodes[index1];
		NetworkNodes n2 = nodes[index2];
		Queue<Integer> findneb = new LinkedBlockingQueue<Integer>();
		boolean[] isused = new boolean[nodes.length];
		for ( Integer i : n1.neberSet){
			findneb.add(i);
			isused[i] = true;
		}
		while( !findneb.isEmpty() ){
			int nebidx = findneb.poll();
			if ( n2.index == nebidx ){
				return true;
			} else {
				for ( Integer j : nodes[nebidx].neberSet ){
					if (!isused[j]){
						findneb.add(j);
						isused[j] = true;
					}
				}
			}
		}
		return false;
		
//		
//		ScaleFree scale=new ScaleFree();
//		return scale.FindItOldParent(nodes, nodes[index1])==scale.FindItOldParent(nodes, nodes[index2]);
	}
//	public static List<NetworkNodes> selectNnodes(NetworkNodes[] nodes,int n){
//		
//	}
	public static Set<Integer> getMostComponentForOneNetwork(NetworkNodes[] nodes)
	{
		int max=0;
		Set<Integer> result = new HashSet<Integer>();
		Set<Integer> component = null;
		for(int i=0;i<nodes.length;i++)
		{
			if ( component == null || !component.contains(i) ){
				component = getComponentByIndex(nodes, nodes[i].index);
			} else {
				continue;
			}
			if ( component.size() > max ){
				max = component.size();
				result = component;
			}
		}
		return result;
	}
	public static int CheckMostComponentForOneNetwork(NetworkNodes[] nodes)
	{
		int max=0;
		boolean[] isTest = new boolean[nodes.length];
		for(int i=0;i<nodes.length;i++)
		{
			if(isTest[i])continue;
//			Boolean[] isConnected = BFS(i,nodes);
			Set<Integer> component = getComponentByIndex(nodes, nodes[i].index);
			int size=0;
			if ( component.size() > max ){
				max = component.size();
			}
			for ( Integer nodeidx : component ){
				isTest[nodeidx] = true;
			}
		}
		return max;
	}
}
