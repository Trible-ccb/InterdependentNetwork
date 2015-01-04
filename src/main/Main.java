package main;
import graphs.ERGraph;
import graphs.Network;
import graphs.NetworkNodes;
import graphs.Price;

import java.io.ByteArrayInputStream;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Map.Entry;




public class Main {
	

	public static String getMemoryInfo() {
	       Runtime currRuntime = Runtime.getRuntime ();
	       int nFreeMemory = ( int ) (currRuntime.freeMemory() / 1024 / 1024);
	       int nTotalMemory = ( int ) (currRuntime.totalMemory() / 1024 / 1024);
	       return nFreeMemory + "M/" + nTotalMemory + "M(free/total)" ;
	}
	
	static int step = 0;//0 代表跑理论值，1代表跑模拟值
	static int endstep = 0;
	static int n = 10000;
	public static void main(String[] args){
		while (step <= endstep ){
//			for ( int i = 4; i <= 4 ; i+=2)
//			runER(n, i);
			for ( float i = 2.8f; i <= 2.8f; i+=0.2f )
			runSF(n, 4,i);
			step+=1;
		}
//		runERRecovery(n, 4);
//		runSFRecovery(n,4, 2.8f);
	}
	public static void runRecovery(Network n1,Network n2,OnGotRecoveryResult mCallback){
		Date start = new Date();
		OnGotRecoveryResult callback = mCallback;
		float initRecovery =  1f;
		float recoveryLinkPro = 0.5f;
		float qc = 0;
		InterdependentOperate.dependentTwoNetRandom(n1.networkNodes, n2.networkNodes);
		for ( initRecovery = 0.1f ; initRecovery <= 0.9f;initRecovery+=0.05f){
			for ( recoveryLinkPro = 0.1f ; recoveryLinkPro <= 0.5f;recoveryLinkPro+=0.05f){
				
				System.out.println("initRecovery=" + initRecovery + " recoveryLinkPro = " +recoveryLinkPro);
				NetworkNodes[] A = copyNetworks(n1.networkNodes);
				NetworkNodes[] B = copyNetworks(n2.networkNodes);
				float cmp = 1f;
				InterdependentOperate.cascadeRecovery(A, B, initRecovery, recoveryLinkPro);
				cmp = (InterdependentOperate.getMostComponentForOneNetwork(A).size()) 
						/ ((float)n1.nodes);
//				if ( cmp < 0.05f && qc == 0 ){
//					qc = p;
//				}
				if ( callback != null ){
					callback.onOneDone(initRecovery, recoveryLinkPro, cmp, n1, n2);
				}
			}
		}
		if ( callback != null ){
			callback.onAllDone(n1, n2,false);
		}
		Date end = new Date();
		System.out.println("operation done.\n"+start.toString()+"\n"+end.toString());
	}
	public static void run(Network n1,Network n2,OnGotOneResult mCallback){
		Date start = new Date();
		OnGotOneResult callback = mCallback;
		float alpha = 1f;
		float beta = 0.5f;
		float fraction = 0.3f;
		
		String filepre1 = n1.displayname; 
		String filepre2 =  n2.displayname;
		
		saveNetDistribution(
				n1,
				filepre1+"Distribution.txt");
		saveNetDistribution(
				n2,
				filepre2+ "Distribution.txt");
		InterdependentOperate.dependentTwoNetRandom(n1.networkNodes, n2.networkNodes);
		int q = 3;
		float delta = 1f;
		for(alpha=-q;alpha<=q;alpha+=delta)
		{
			for(beta=-q;beta<=q;beta+=delta)
			{
				float pc = 0;
				HashMap<Integer, Double> faiKs = null;
				HashMap<Integer, Double> faiKs2 = null;
				if( step == 0 ){
					faiKs = InterdependentOperate.funFAI(n1, n2, alpha, beta);
					faiKs2 = InterdependentOperate.funFAI(n2,n1,beta, alpha);
				}
				for(float i=0.1f;i<=0.9f;i+=0.05f)
				{
					
					fraction = i;
					float p = 1 - fraction;
					boolean isTheory;
					float cmp = 0;
					System.out.println("alpha=" + alpha + " beta=" + beta + " 1-p=" + fraction);
					if ( step == 0 ){
//						计算理论值
						isTheory = true;
						double tcmp = InterdependentOperate.cascadeDependentNetInTheory(
								faiKs,faiKs2,n1, n2, fraction);
						cmp = (float) tcmp;
					} else {
						isTheory = false;
						NetworkNodes[] A = copyNetworks(n1.networkNodes);
						NetworkNodes[] B = copyNetworks(n2.networkNodes);
						InterdependentOperate.cascadeDependentNet(A, B, alpha, fraction,beta);
						cmp = (InterdependentOperate.CheckMostComponentForOneNetwork(A)) 
								/ ((float)n1.nodes*p);
					}
					if ( cmp < 0.1f && pc == 0 ){
						pc = p;
					}
					if ( callback != null ){
						callback.onOneDone(alpha, beta, p,pc, cmp, n1, n2,isTheory);
					}
				}
			}
		}
		if ( callback != null ){
			callback.onAllDone(n1, n2);
		}
		Date end = new Date();
		System.out.println("operation done.\n"+start.toString()+"\n"+end.toString());
	}
	public static void runSF(int n,int avgdeg,final float lamda){
		int num = n;
		Price price1 = new Price();
		Network sf1 = null;
		Network sf2 = null;
		int m = avgdeg / 2;
		float a = (lamda - 2) * m;
		sf1 = price1.buildNetwork(num, avgdeg / 2, a);
		Price price2 = new Price();
		sf2 = price2.buildNetwork(num, avgdeg / 2, a);
		OnGotOneResult callback = new OnGotOneResultImpl(lamda);
		OnGotOneResultImpl.avgKVsPc = new LinkedHashMap<Object, Map<Object,Object>>();
		run(sf1, sf2,callback);
	}
	public static void runSFRecovery(int n,int avgdeg,final float lamda){
		int num = n;
		Price price1 = new Price();
		Network sf1 = null;
		Network sf2 = null;
		int m = avgdeg / 2;
		float a = (lamda - 2) * m;
		sf1 = price1.buildNetwork(num, avgdeg / 2, a);
		Price price2 = new Price();
		sf2 = price2.buildNetwork(num, avgdeg / 2, a);
		OnGotRecoveryResult callback = new OnGotRecoveryResultImpl(lamda);
//		OnGotOneResultImpl.avgKVsPc = new LinkedHashMap<Object, Map<Object,Object>>();
		runRecovery(sf1, sf2,callback);
	}
	public static void runSFAndER(int n,int avgdeg,float lamda){
		int num = n;
		Price price1 = new Price();
		ERGraph er = new ERGraph();
		Network n1 = null;
		Network n2 = null;
		int m = avgdeg / 2;
		float a = (lamda - 2) * m;
		n1 = price1.buildNetwork(num, avgdeg / 2, a);
		n2 = er.buildERGraph1(num, avgdeg);
		OnGotOneResult callback = new OnGotOneResultImpl(lamda);
		OnGotOneResultImpl.avgKVsPc = new LinkedHashMap<Object, Map<Object,Object>>();
		run(n1, n2,callback);
	}
	public static void runERAndSF(int n,int avgdeg,float lamda){
		int num = n;
		Price price1 = new Price();
		ERGraph er = new ERGraph();
		Network n1 = null;
		Network n2 = null;
		int m = avgdeg / 2;
		float a = (lamda - 2) * m;
		n1 = price1.buildNetwork(num, avgdeg / 2, a);
		n2 = er.buildERGraph1(num, avgdeg);
		OnGotOneResult callback = new OnGotOneResultImpl(lamda);
		OnGotOneResultImpl.avgKVsPc = new LinkedHashMap<Object, Map<Object,Object>>();
		run(n2, n1,callback);
	}
	public static void runER(int n,int avgdeg){
		Network net1 = null;
		Network net2 = null;
		ERGraph er = new ERGraph();
		net1 = er.buildERGraph1(n, avgdeg);
		ERGraph er2 = new ERGraph();
		net2 = er2.buildERGraph1(n, avgdeg);
		OnGotOneResult callback = new OnGotOneResultImpl(0);
		if(OnGotOneResultImpl.avgKVsPc == null){
			OnGotOneResultImpl.avgKVsPc = new LinkedHashMap<Object, Map<Object,Object>>();
		}
		run(net1, net2,callback);
		
	}
	public static void runERRecovery(int n,int avgdeg){
		Network net1 = null;
		Network net2 = null;
		ERGraph er = new ERGraph();
		net1 = er.buildERGraph1(n, avgdeg);
		ERGraph er2 = new ERGraph();
		net2 = er2.buildERGraph1(n, avgdeg);
		OnGotRecoveryResult callback = new OnGotRecoveryResultImpl(0);
//		if(OnGotOneResultImpl.avgKVsPc == null){
//			OnGotOneResultImpl.avgKVsPc = new LinkedHashMap<Object, Map<Object,Object>>();
//		}
		runRecovery(net1, net2,callback);
	}
	public static void saveResults(String filename,String str){
		StringBuffer sb = new StringBuffer();
//		for ( int i = 0 ; i < nodes.length ; i++ ){
//			sb.append(nodes[i].toString());
//		}
		FileUtils.storeInputStream(
				filename, new ByteArrayInputStream(str.getBytes()), "utf-8");
	}
	
	public static void saveNetDistribution(Network net,String filename){
		System.out.println("saving Distribution:" + filename);
		StringBuffer sb = new StringBuffer();
		HashMap<Integer, Float> degrees = net.getDegreeDistribution();
		for ( Entry<Integer, Float> e : degrees.entrySet() ){
			sb.append(e.getKey() + " " + e.getValue() + "\n");
		}
		FileUtils.storeInputStream(
				filename, new ByteArrayInputStream(sb.toString().getBytes()), "utf-8");
		System.out.println("saving Distribution:" + filename + " end");
	}
	public static NetworkNodes[] copyNetworks(NetworkNodes[] net1){
		NetworkNodes[] A = new NetworkNodes[net1.length];
		for(int n1=0;n1<net1.length;n1++)
		{
			A[n1]=new NetworkNodes();
			A[n1].copyValueBy(net1[n1]);
		}
		return A;
	}
	
	static class OnGotOneResultImpl implements OnGotOneResult{
		
		float lamda = 0;
		String preFilenameOfResult = "";//如果是理论值，加个T在文件名头
		
		public OnGotOneResultImpl(float lamdaForSF) {
			lamda = lamdaForSF;
			if ( lamdaVsPc == null )
			lamdaVsPc = new LinkedHashMap<Object, Map<Object, Object>>();
		}
//		保存数据给画图
//		p值为x轴，cmp为y轴
//		outer key = "p" value = Map : inner key "a_b",value cmp
		Map<Object, Map<Object, Object>> pVsCmp = new LinkedHashMap<Object, Map<Object, Object>>();
		
//		lamda值为x轴，pc为y轴
//		outer key = "lamda" value = Map : inner key "a_b",value pc
		static Map<Object, Map<Object, Object>> lamdaVsPc = null;
		
//		<k>值为x轴，pc为y轴
//		outer key = "<k>" value = Map : inner key "a_b",value pc
		public static Map<Object, Map<Object, Object>> avgKVsPc = null;
		
//		a值为x轴，cmp为y轴
//		outer key:alpha value Map: inner key : beta,value cmp
		Map<Object, Map<Object, Object>> abVsPc = new LinkedHashMap<Object, Map<Object,Object>>();

		@Override
		public void onOneDone(float a, float b, float p, float pc, float cmp,
				Network n1, Network n2,boolean isTheory) {
			if( isTheory ){
				preFilenameOfResult = "T";
			}
			Map<Object, Object> abcmp = pVsCmp.get(p);
			if (abcmp == null){
				abcmp = new LinkedHashMap<Object, Object>();
			}
			pVsCmp.put(p, abcmp);
			abcmp.put(a+"_"+b, cmp);
			
			Map<Object, Object> bcmp = abVsPc.get(a);
			if ( bcmp == null ){
				bcmp = new LinkedHashMap<Object, Object>();
			}
			bcmp.put(b, pc);
			abVsPc.put(a, bcmp);
			
			if ( lamda != 0 ){
				Map<Object, Object> abpc = lamdaVsPc.get(lamda);
				if (abpc == null){
					abpc = new LinkedHashMap<Object, Object>();
				}
				lamdaVsPc.put(lamda, abpc);
				abpc.put(a+"_"+b, pc);
			}
			
			Map<Object, Object> abpc = avgKVsPc.get(n1.avgdegree);
			if (abpc == null){
				abpc = new LinkedHashMap<Object, Object>();
			}
			avgKVsPc.put(n1.avgdegree, abpc);
			abpc.put(a+"_"+b, pc);
		}
		private void resolveMap(String prefirstline,Map<Object,Map<Object, Object>> map,String filename){
			StringBuffer firstline = new StringBuffer().append(prefirstline);
			StringBuffer valueline = new StringBuffer();
			boolean readfirst = false;
			for ( Entry<Object,Map<Object, Object>> entry : map.entrySet() ){
				Map<Object, Object> value = entry.getValue();
				Object keyp = entry.getKey();
				valueline.append("\n").append(keyp);
				for ( Entry<Object, Object> entry2 : value.entrySet() ){
					Object keyab = entry2.getKey();
					Object cmp = entry2.getValue();
					if(!readfirst)
					firstline.append("\t").append(keyab);
					valueline.append("\t").append(cmp);
				}
				readfirst = true;
			}
			saveResults(preFilenameOfResult+filename, (firstline.append(valueline)).toString());
		}
		@Override
		public void onAllDone(Network n1, Network n2) {
			String filepre1 = n1.displayname; 
			String filepre2 =  n2.displayname;
			String filenamepVsCmp = filepre1+filepre2+"_pVsCmp.txt";
			String filenameLamdaVsPc = filepre1+filepre2+"_LamdaVsPc.txt";
			String filenameKVsPc = filepre1+filepre2+"_avgKVsPc.txt";
			String filenameabVsPc = filepre1+filepre2+"_abVsPc.txt";
			
			resolveMap("p", pVsCmp, filenamepVsCmp);
			resolveMap("a", abVsPc, filenameabVsPc);
			resolveMap("<k>", avgKVsPc, filenameKVsPc);
			if ( !lamdaVsPc.isEmpty() ){
				resolveMap("lamda", lamdaVsPc, filenameLamdaVsPc);
			}
		}
	}
	static class OnGotRecoveryResultImpl implements OnGotRecoveryResult{
		float lamda = 0;
		String preFilenameOfResult = "Rec";
		public OnGotRecoveryResultImpl(float lamda){
			this.lamda = lamda;
		}
//		保存数据给画图
//		p值为x轴，cmp为y轴
//		outer key = "p" value = Map : inner key "a_b",value cmp
		Map<Object, Map<Object, Object>> qVsCmp = new LinkedHashMap<Object, Map<Object, Object>>();
		@Override
		public void onOneDone(float initRecovery, float recoveryLinkPro,
				float cmp, Network n1, Network n2) {
			Map<Object, Object> rcmp = qVsCmp.get(initRecovery);
			if (rcmp == null){
				rcmp = new LinkedHashMap<Object, Object>();
			}
			qVsCmp.put(initRecovery, rcmp);
			rcmp.put(recoveryLinkPro+"", cmp);
		}

		@Override
		public void onAllDone(Network n1, Network n2,boolean isTheory) {
			if( isTheory ){
				preFilenameOfResult += "T";
			}
			String filepre1 = n1.displayname; 
			String filepre2 =  n2.displayname;
			String filenamepVsCmp = filepre1+filepre2+"_qVsCmp.txt";
			resolveMap("q", qVsCmp, filenamepVsCmp);
		}
		private void resolveMap(String prefirstline,Map<Object,Map<Object, Object>> map,String filename){
			StringBuffer firstline = new StringBuffer().append(prefirstline);
			StringBuffer valueline = new StringBuffer();
			boolean readfirst = false;
			for ( Entry<Object,Map<Object, Object>> entry : map.entrySet() ){
				Map<Object, Object> value = entry.getValue();
				Object keyp = entry.getKey();
				valueline.append("\n").append(keyp);
				for ( Entry<Object, Object> entry2 : value.entrySet() ){
					Object keyab = entry2.getKey();
					Object cmp = entry2.getValue();
					if(!readfirst)
					firstline.append("\t").append(keyab);
					valueline.append("\t").append(cmp);
				}
				readfirst = true;
			}
			saveResults(filename, (firstline.append(valueline)).toString());
		}
	}
	public interface OnGotOneResult{
		
		/**
		 * @param a alpha
		 * @param b beta
		 * @param p = 1 - fraction 
		 * @param pc  
		 * @param cmp 序参量
		 */
		void onOneDone(float a,float b,float p,float pc,float cmp,Network n1,Network n2,boolean isTheory);
		
		void onAllDone(Network n1,Network n2);
	}
	public interface OnGotRecoveryResult{
		void onOneDone(float initRecovery,float recoveryLinkPro,float cmp,Network n1,Network n2);
		void onAllDone(Network n1,Network n2,boolean isTheory);
	}
}
