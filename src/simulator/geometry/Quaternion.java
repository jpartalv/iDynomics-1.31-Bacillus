package simulator.geometry;

/** 
 * @author  Jean Pierre
 */
public class Quaternion
{
    public double w = 0;
    public double x = 0;
    public double y = 0;
    public double z = 0;
    
    public Quaternion(double w, double x, double y, double z)
    {
        this.w = w;
        this.x = x;
        this.y = y;
        this.z = z;
    }

    /*public Quaternion(float w, Vector v)
    {
        this.w = w;
        this.x = Convert.ToSingle(v.x);
        this.y = Convert.ToSingle(v.y);
        this.z = Convert.ToSingle(v.z);
    }*/

    public static Quaternion Multiply(Quaternion a, Quaternion b)
    {
    	double nw = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z;
    	double nx = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y;
    	double ny = a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z;
    	double nz = a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x;
        return new Quaternion(nw, nx, ny, nz);

    }

    public static Quaternion Conjugate(Quaternion a)
    {
        return new Quaternion(a.w, -a.x, -a.y, -a.z);
    }
}