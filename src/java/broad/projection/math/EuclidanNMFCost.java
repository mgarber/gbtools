package broad.projection.math;


public class EuclidanNMFCost implements NMFCostFunction {

	public double cost(Matrix original, Matrix approximation) {
		return original.euclidean(approximation);
	}
	
	public Matrix[] update(Matrix V, Matrix W, Matrix H){
		Matrix Hnew = updateH(V,W,H);
		Matrix Wnew = updateW(V,W,Hnew);
		Matrix[] rtrn={Wnew, Hnew};
		return rtrn;
	}
	
	public Matrix updateH(Matrix V, Matrix W, Matrix H){
		Matrix WT=W.ntranspose();
		Matrix WTW=WT.nmatmul(W);
		Matrix denom= WTW.nmatmul(H);
		Matrix num=WT.nmatmul(V);
		Matrix quot=num.ndiv(denom);
		Matrix Hnew=H.nmul(quot);
		return Hnew;
	}
	
	public Matrix updateW(Matrix V, Matrix W, Matrix H){
		Matrix HT=H.ntranspose();
		Matrix num=V.nmatmul(HT);
		Matrix WH=W.nmatmul(H);
		Matrix denom=WH.nmatmul(HT);
		Matrix quot=num.ndiv(denom);
		Matrix Wnew=W.nmul(quot);
		return Wnew;
	}
	
	public String getName() {return "Euclidean";}
}
