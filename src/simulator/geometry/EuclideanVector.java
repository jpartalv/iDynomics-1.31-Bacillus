package simulator.geometry;
//import java.awt.Point;

/** 
 * @author  Jean Pierre
 */
public class EuclideanVector
{
	public double[] start =  {0,0,0};
    public double[] end = {0,0,0};
    public double magnitude;
    public double mag_x;
    public double mag_y;
    public double mag_z;

   
    
    public EuclideanVector(ContinuousVector start, ContinuousVector end)
    {
    	this.start[0] = start.x;
    	this.start[1] = start.y;
    	this.start[2] = start.z;
        this.end[0] = end.x;
        this.end[1] = end.y;
        this.end[2] = end.z;
        mag_x = end.x - start.x;
        mag_y = end.y - start.y;
        mag_z = end.z - start.z;
        
        magnitude = (float)(Math.sqrt(Math.pow(mag_x, 2) + Math.pow(mag_y, 2) + Math.pow(mag_z, 2)));
    }
    
    public EuclideanVector(double[] start, double[] end)
    {
    	this.start = start;
        this.end = end;
        mag_x = end[0] - start[0];
        mag_y = end[1] - start[1];
        mag_z = end[2] - start[2];
        
        magnitude = (float)(Math.sqrt(Math.pow(mag_x, 2) + Math.pow(mag_y, 2) +  + Math.pow(mag_z, 2)));
    }

        public EuclideanVector(double[] start, double xMag, double yMag, double zMag)
        {
            this.start = start;
            this.end[0] = start[0] + xMag;
            this.end[1] = start[1] + yMag;
            this.end[2] = start[2] + zMag;
            
            mag_x = xMag;
            mag_y = yMag;
            mag_z = zMag;
            magnitude = Math.sqrt(Math.pow(mag_x, 2) + Math.pow(mag_y, 2) +  + Math.pow(mag_z, 2));
        }

    
        /*public Vector2D Clone()
        {
           Vector clone = new Vector();
           clone.start = this.start;
           clone.end = this.end;
           clone.magnitude = this.magnitude;
           clone.mag_x = this.mag_x;
           clone.mag_y = this.mag_y;
           clone.mag_z = this.mag_z;
           return clone;
        }*/

        /*protected void CalcMagnitudes()
        {
            this.mag_x = endPoint.Xf - startPoint.Xf;
            this.mag_y = endPoint.Yf - startPoint.Yf;
            this.mag_z = endPoint.Zf - startPoint.Zf;
            magnitude = Convert.ToSingle(Math.Sqrt(Math.Pow(mag_x, 2) + Math.Pow(mag_y, 2) + Math.Pow(mag_z, 2)));
        }*/

        /*public void Reverse()
        {
            Point3 temp  = this.start;
            this.start = this.endPoint;
            this.end = temp;
            this.CalcMagnitudes();
        }*/

        /* operations */
        public EuclideanVector Plus(EuclideanVector B)
        {
            mag_x += B.mag_x;
            mag_y += B.mag_y;
            mag_z += B.mag_z;
            return new EuclideanVector(start,mag_x,mag_y,mag_z);
            /*this.end = new Point();
            end.setLocation(start.getX() + mag_x, start.getY()+ mag_y);
            magnitude = (float)(Math.sqrt(Math.pow(mag_x, 2) + Math.pow(mag_y, 2) + Math.pow(mag_z, 2)));*/
        }
        
        public EuclideanVector Plus(double x, double y, double z)
        {
            mag_x += x;
            mag_y += y;
            mag_z += z;
            return new EuclideanVector(start,mag_x,mag_y,mag_z);
            /*this.end = new Point();
            end.setLocation(start.getX() + mag_x, start.getY()+ mag_y);
            magnitude = (float)(Math.sqrt(Math.pow(mag_x, 2) + Math.pow(mag_y, 2) + Math.pow(mag_z, 2)));*/
        }    
        
        public EuclideanVector Minus(EuclideanVector B)
        {
            mag_x -= B.mag_x;
            mag_y -= B.mag_y;
            mag_z -= B.mag_z;
            return new EuclideanVector(start,mag_x,mag_y,mag_z);
            /*this.end = new Point();
            end.setLocation(start.getX() + mag_x, start.getY()+ mag_y);
            magnitude = (float)(Math.sqrt(Math.pow(mag_x, 2) + Math.pow(mag_y, 2) + Math.pow(mag_z, 2)));*/
        }
        
   	 	public EuclideanVector Times(double a) // escalar and vector
   	 	{
   	 		return new EuclideanVector(start, mag_x * a, mag_y * a, mag_z * a);
   	 	}

        public double DotProduct(EuclideanVector B)
        {
            return ((mag_x * B.mag_x) + (mag_y * B.mag_y) + (mag_z * B.mag_z));
        } 
        
        public double DotProduct(ContinuousVector B)
        {
            return ((mag_x * B.x) + (mag_y * B.y) + (mag_z * B.z));
        } 
        
        public EuclideanVector CrossProduct(EuclideanVector B)
        {
            
        	EuclideanVector v = new EuclideanVector(start,((mag_y * B.mag_z) - (mag_z * B.mag_y)),
                    ((mag_z * B.mag_x) - (mag_x * B.mag_z)),
                    ((mag_x * B.mag_y) - (mag_y * B.mag_x)));
            return v;
        }
        
        public EuclideanVector Normalize()
        {
        	return new EuclideanVector(start,mag_x / magnitude,mag_y / magnitude,mag_z / magnitude);
        }
        
        public static EuclideanVector rescaleVector(ContinuousVector[] points, double newLength)
        {
        	double magnitude = points[1].distance(points[0]);
        	double magX = ((points[1].x - points[0].x) / magnitude) * newLength;
			double magY = ((points[1].y - points[0].y) / magnitude) * newLength;
			double magZ = ((points[1].z - points[0].z) / magnitude) * newLength;
			points[1].x = points[0].x + magX;
			points[1].y = points[0].y + magY;
			points[1].z = points[0].z + magZ;
        	
			return new EuclideanVector(points[0],points[1]);
        }
        
        public ContinuousVector getContinuousVector()
        {
        	return new ContinuousVector(mag_x,mag_y,mag_z);
        }
    }


