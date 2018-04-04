package broad.projection.math;

public interface NMFCostFunction {
	public double cost(Matrix original, Matrix approximation);
	public String getName();
	public Matrix[] update(Matrix V, Matrix W, Matrix H);
	//public Matrix updateW(Matrix V, Matrix W, Matrix H);
	//public Matrix updateH(Matrix V, Matrix W, Matrix H);
}
