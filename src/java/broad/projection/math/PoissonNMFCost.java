package broad.projection.math;

public class PoissonNMFCost implements NMFCostFunction {

	public double cost(Matrix original, Matrix approximation) {
		return original.probabilistic(approximation);
	}
	public String getName() { return "Poisson";}

	
	public Matrix[] update(Matrix V, Matrix W, Matrix H){
		Matrix WH=W.nmatmul(H); 
		Matrix WT=W.ntranspose();
		Matrix K = V.ndiv(WH); //interm
		Matrix HT=H.ntranspose();
		
		//HNew
		Matrix num = WT.nmatmul(K);
		Matrix Hnew = H.nmul(num);
		
		Matrix sum = W.colSum();  // Sum the columns of W

        // LOOP THROUGH EACH ROW IN H
        for (int i = 0; i < H.nr; i = i + 1) {
            // LOOP THROUGH EACH COLUMN IN H
            for (int j = 0; j < H.nc; j = j + 1) {
                // NORMALIZE THE ROW
                Hnew.x[i][j] = (Hnew.x[i][j]) / (sum.x[0][i]);
            }
        }
        
        H=Hnew;
        
        //WNew
        Matrix Wnum = K.nmatmul(HT);
		Matrix Wnew = W.nmul(Wnum);
		
		Matrix Wsum = H.rowSum();  // Sum the columns of H
		
		 // LOOP THROUGH EACH COLUMN IN W
        for (int j = 0; j < W.nc; j = j + 1) {
            // LOOP THROUGH EACH ROW IN W
            for (int i = 0; i < W.nr; i = i + 1) {
                // NORMALIZE THE ROW
                Wnew.x[i][j] = (Wnew.x[i][j]) / (Wsum.x[0][j]);
            }
        }
        
        Matrix[] rtrn={Wnew, Hnew};
        return rtrn;
	}
	
	public Matrix updateH(Matrix V, Matrix W, Matrix H){
		Matrix WH=W.nmatmul(H); 
		Matrix WT=W.ntranspose();
		Matrix K = V.ndiv(WH); //interm
		
		Matrix num = WT.nmatmul(K);
		Matrix Hnew = H.nmul(num);
		
		Matrix sum = W.colSum();  // Sum the columns of W

        // LOOP THROUGH EACH ROW IN H
        for (int i = 0; i < H.nr; i = i + 1) {
            // LOOP THROUGH EACH COLUMN IN H
            for (int j = 0; j < H.nc; j = j + 1) {
                // NORMALIZE THE ROW
                Hnew.x[i][j] = (Hnew.x[i][j]) / (sum.x[0][i]);
            }
        }
        return Hnew;
	}
	
	public Matrix updateW(Matrix V, Matrix W, Matrix H){
		Matrix WH=W.nmatmul(H); 
		Matrix HT=H.ntranspose();
		Matrix K = V.ndiv(WH);
		Matrix num = K.nmatmul(HT);
		Matrix Wnew = W.nmul(num);
		
		Matrix sum = H.rowSum();  // Sum the columns of H
		
		 // LOOP THROUGH EACH COLUMN IN W
        for (int j = 0; j < W.nc; j = j + 1) {
            // LOOP THROUGH EACH ROW IN W
            for (int i = 0; i < W.nr; i = i + 1) {
                // NORMALIZE THE ROW
                Wnew.x[i][j] = (Wnew.x[i][j]) / (sum.x[0][j]);
            }
        }
		return Wnew;
	}
	
}
