public class complex 
{
    private double a; // real value for objects of type complex
    private double b; // imaginary value for objects of type complex
    

    /*
     * constructor
     * constructs objects of type complex
     * a and b are respectively real and imaginary values of complex number
     * four constructor methods
     * 2 double, 1 double, 1 complex, and no parameter methods
     * return void: creates complex objects
     */
    public complex(double x, double y)
    {
        a = x;
        b = y;
    }
    public complex(double x){
        this(x, 0);
    }
    public complex(complex x){
        a = x.a;
        b = x.b;
    }
    public complex(){
        this(0, 0); // complex with real and imaginary values (0,0)
    }

    
    /*
     * getters
     * gets real and imaginary values of complex objects
     * two nonstatic methods
     * returns real and imaginary values for given complex value as type double
     */
    public double getReal(){
        return(a);
    }
    public double getImaginary(){
        return(b);
    }


    /*
     * setters
     * sets values for real and imaginary values of complex objects
     * two nonstatic methods
     * parameter: double val set as real or imaginary
     * return void: changes values of a complex object
     */
    public void setReal(double val){
        a = val;
    }
    public void setImaginary(double val){
        b = val;
    }


    /*
     * arg
     * finds arctan in radians of a complex number in range (0, 2 * PI)
     * static and nonstatic methods available
     * static parameter complex x is the complex by which arctan is taken from
     * returns sol: new object of type double that holds arctan of complex value
     */
    public double arg(){
        if(a == 0 & b == 0){
            return 9999999;
        }
        double sol = Math.atan2(b, a); // stores arctan of complex value in sol (short for solution)
        sol += (Math.PI * 2);
        sol = sol % (Math.PI * 2); // adjusts range of arctan from (-PI, PI) to (0, 2 * PI)
        return sol;
    }
    public static double arg(complex x){
        return x.arg();
    }


    /*
     * mag
     * finds distance from 0 of complex value
     * static and nonstatic methods available
     * static parameter complex x is the complex that's distance from 0 is found
     * returns sol: new object of type double that holds value of distance from 0 in type double
     */
    public double mag(){
        double sol = Math.sqrt(a * a + b * b); // calculations for distance from 0 (hypotenuse where real is x axis and imaginary is y axis)
        return sol;
    }
    public static double mag(complex x){
        return x.mag();
    }


    /*
     * conjugate
     * finds conjugate of given complex number
     * static and nonstatic methods available
     * static parameter complex x is the complex that's conjugate is found
     * returns sol: new complex object that holds the conjugate
     */
    public complex conjugate(){
        complex sol = new complex(a, -b); // conjugate has negative imaginary value
        return sol;
    }
    public static complex conjugate(complex x){
        return x.conjugate();
    }


    /*
     * add
     * adds two complex values
     * static and nonstatic methods available
     * nonstatic parameter: complex y is value that is added
     * static parameter complex x and complex y are the values that are added
     * returns sol: new complex object that holds the sum of the two complex values
     */
    public complex add(complex y){
        double solA = a + y.a; // adds reals
        double solB = b + y.b; // adds imaginary
        complex sol = new complex(solA, solB); // creates new object with solution
        return sol;
    }
    public static complex add(complex x, complex y){
        return x.add(y);
    }


    /*
     * subtract
     * subtracts two complex values
     * parameter: complex y is value that is subtracted
     * static and nonstatic methods available
     * static parameter complex y is subtracted from parameter complex x
     * returns sol: new complex object that holds the difference of the two complex values
     */
    public complex subtract(complex y){
        double solA = a - y.a; // subtracts reals
        double solB = b - y.b; // subtracts imaginary
        complex sol = new complex(solA, solB);
        return sol;
    }
    public static complex subtract(complex x, complex y){
        return x.subtract(y);
    }
    

    /*
     * divide
     * divides two complex values
     * static and nonstatic methods available
     * nonstatic parameter: complex y is the denominator
     * static parameter: complex x and complex y are the numerator and denominator respectively
     * returns sol: new complex object that holds the quotient of the two complex values
     */
    public complex divide(complex y){
        if(y == new complex()){
            return new complex(9999999);
        }
        double den = Math.pow(y.a, 2) + Math.pow(y.b, 2); // defines the denominator
        double solA = ((a * (y.a)) - (b * y.conjugate().b)) / den; // defines solA (real part of the final solution)
        double solB = ((a * y.conjugate().b) + (b * y.a)) / den; // defines solB (imaginary part of the final solution)
        complex sol = new complex(solA, solB);
        return sol;
    }
    public static complex divide(complex x, complex y){
        return x.divide(y);
    }


    /*
     * multiply
     * multiplies two complex values
     * static and nonstatic methods available
     * nonstatic parameter: complex y is value that is multiplied
     * static parameter: complex x and complex y are the values that are multiplied
     * returns sol: new complex object that holds the product of the two complex values
     */
    public complex multiply(complex y){
        double solA = ((a * (y.a)) - (b * y.b)); // calculations for solA
        double solB = ((a * y.b) + (b * y.a)); // calculations for solB
        complex sol = new complex(solA, solB);
        return sol;
    }
    public static complex multiply(complex x, complex y){
        return x.multiply(y);
    }


    /*
     * toString
     * converts complex value into readable string format
     * returns sol: new string object with readable format of complex value
     */
    public String toString(){
        String sol = "(" + a + ", " + b + "i)"; // conversion to string
        return sol;
    }


    /*
     * equals
     * determines if two complex values are equivlant
     * static and nonstatic methods available
     * nonstatic parameter: complex y is value that is compared
     * static parameter: complex x and complex y are the values that are compared
     * returns boolean: boolean object indicating whether objects are equivlant (true) or not (false)
     */
    public boolean equals(complex y, double thresh){
        double epsilonA; // epsilon of real part
        double epsilonB; // epsilon of imaginary part
        if(a != 0){
            epsilonA = (a - y.a) / a;
        }
        else if(y.a != 0){
            epsilonA = (a - y.a) / y.a;
        }
        else{
            epsilonA = 0;
        }
        if(b != 0){
            epsilonB = (b - y.b) / b;
        }
        else if(y.b != 0){
            epsilonB = (b - y.b) / y.b;
        }
        else{
            epsilonB = 0;
        }
        if(epsilonA < thresh && epsilonA > -thresh && epsilonB < thresh && epsilonB > -thresh){ // condition for whether the values are equivlant
            return true;
        }
        else if((y.a-a) <= thresh && (y.a-a) >= -thresh && (y.b-b) <= thresh && (y.b-b) >= -thresh && (a==0 || y.a==0 || b==0 || y.b==0)){
            return true;
        }
        else{
            return false;
        }
    }
    public boolean equals(complex y){
        return equals(y, 0.000001);
    }
    public static boolean equals(complex x, complex y, double thresh){
        return x.equals(y, thresh);
    }
    public static boolean equals(complex x, complex y){
        return x.equals(y, 0.000001);
    }


    /*
     * complexequals
     * compares complex numbers using equals method
     * static and nonstatic methods available
     * nonstatic parameter: object o is value that is compared
     * static parameter: complex x and object o are the values that are compared
     * returns boolean: boolean object indicating whether objects are equivlant (true) or not (false)
     */
    public boolean complexequals(Object o, double thresh){
        if (this.getClass() == o.getClass()){
            complex w = (complex) o;
            return equals(w, thresh);
        }
        else{
            return false;
        }
    }
    public boolean complexequals(Object o){
        if (this.getClass() == o.getClass()){
            complex w = (complex) o;
            return equals(w);
        }
        else{
            return false;
        }
    }
    public static boolean complexequals(complex x, Object o, double thresh){
        return x.complexequals(o, thresh);
    }
    public static boolean complexequals(complex x, Object o){
        return x.complexequals(o);
    }


    /*
     * power
     * raises a complex value to a double power
     * static and nonstatic methods available
     * nonstatic parameter: double y is the value to which the complex is raised
     * static parameter: complex x is the complex that will be raised to the power double y
     * returns sol: new complex object that holds complex value raised to power double y
     */
    public complex power(double y){
        double M = mag();
        double theta = arg();
        double solA = Math.pow(M, y) * Math.cos(y * theta); // calculations for solA
        double solB = Math.pow(M, y) * Math.sin(y * theta); // calculations for solB
        complex sol = new complex(solA, solB);
        return sol;
    }
    public static complex power(complex x, double y){
        return x.power(y);
    }
    

    /*
     * kroot
     * finds all n possible nth roots of a given complex value
     * static and nonstatic methods available
     * nonstatic parameter: int n represents which root
     * static parameter: complex x is the complex that's root is taken and int n represents which root
     * returns arr_sol: new complex array which holds all n nth roots
     */
    public complex[] kroot(int n){ 
        complex[] arr_sol = new complex[n]; // n possible solutions for nth root in arr_sol (short for "array of solutions")
        double R = mag();
        double theta = arg();
        double solA;
        double solB;
        double N = (double) n;
        for(int k = 0; k < n; k++){
            double K = (double) k;
            solA = Math.pow(R, (1/N));
            solA *= Math.cos((theta + (2 * Math.PI * K)) / N); // calculations for solA (real part)
            solB = Math.pow(R, (1/N));
            solB *= Math.sin((theta + (2 * Math.PI * K)) / N); // calculations for isolB (imaginary part)
            arr_sol[k] = new complex(solA, solB); // real and imaginary value stored as complex value in complex array arr_sol[]
        }
        return arr_sol;
    }
    public static complex[] kroot(complex x, int n){
        return x.kroot(n);
    }


    /*
     * complexpower
     * raises a complex value to a complex power
     * static and nonstatic methods available
     * nonstatic parameter: complex y is the value to which the complex is raised
     * static parameter: complex x is the complex that will be raised to the power complex y
     * returns sol: new complex object that holds complex value raised to power complex y
     */
    public complex complexpower(complex y){
        double r = mag();
        double theta = arg();
        complex gB = new complex(Math.log(r), theta);
        complex gA = y;
        complex g = gB.multiply(gA); // solve for complex g used in main calculations to solve for solA and solB
        double solA = Math.exp(g.a) * Math.cos(g.b); // calculations for solA
        double solB = Math.exp(g.a) * Math.sin(g.b); // calculations for solB
        complex sol = new complex(solA, solB);
        return sol;
    }
    public static complex complexpower(complex x, complex y){
        return x.complexpower(y);
    }


    /*
     * complexln
     * finds the ln of a complex number
     * static and nonstatic methods available
     * static parameter: complex x, ln of complex x will be taken
     * returns sol: new complex object that holds ln of complex value
     */
    public complex complexln(){
        double solA = Math.log(mag()); // calculations for solA 
        double solB = arg(); // calculations for solB
        complex sol = new complex(solA, solB);
        return sol;
    }
    public static complex complexln(complex x){
        return x.complexln();
    }


    /*
     * complexsin
     * finds the sine of a complex value
     * static and nonstatic methods available
     * static parameter: complex x is the complex that's sine will be found
     * returns sol: new complex object that holds sine of complex value
     */
    public complex complexsin(){
        double solA1 = Math.exp(-b) * Math.sin(a); // double solA1 and solA2 used to find solA
        double solB1 = Math.exp(-b) * Math.cos(a); // double solB1 and solB2 used to find solB
        double solA2 = Math.exp(b) * Math.sin(-a);
        double solB2 = Math.exp(b) * Math.cos(a);
        double solA = 0.5 * (-solA1 + solA2); // calculations for solA
        double solB = 0.5 * (-solB1 + solB2); // calculations for solB
        complex sol = new complex(solA, solB);
        return sol;
    }
    public static complex complexsin(complex x){
        return x.complexsin();
    }


    /*
     * complexcos
     * finds the cosine of a complex number
     * static and nonstatic methods available
     * static parameter: complex x is the complex that's cosine will be found
     * returns sol: new complex object that holds cosine of complex value
     */
    public complex complexcos(){
        double solA1 = Math.exp(-b) * Math.cos(a); // double solA1 and solA2 used to find solA (real part)
        double solB1 = Math.exp(-b) * Math.sin(a); // double solB1 and solB2 used to find solB (imaginary part)
        double solA2 = Math.exp(b) * Math.cos(a);
        double solB2 = Math.exp(b) * Math.sin(a);
        double solA = 0.5 * (solA1 + solA2); // calciulations for solA
        double solB = 0.5 * (solB1 - solB2); // calculations for solB
        complex sol = new complex(solA, solB);
        return sol;
    }
    public static complex complexcos(complex x){
        return x.complexcos();
    }


    /*
     * complexarcsin
     * finds the arcsine of a given complex value
     * static and nonstatic methods available
     * static parameter: complex x is the complex that's arcsin will be found
     * returns arr_sol: new complex array object that holds the two solutions of arcsine of a complex value
     */
    public complex[] complexarcsin(){
        complex xsquare = power(2); 
        complex solA = new complex(1 - xsquare.a, -xsquare.b);
        complex solB = new complex(-b, a);
        complex[] rootA = solA.kroot(2);
        complex[] arr_sol = new complex[2]; // creating array ar_sol (short for "array of solutions") of length 2
        complex sumAB = new complex();
        complex lnSum = new complex();
        for(int i = 0; i < rootA.length; i++){
            sumAB = solB.add(rootA[i]); // calculations for complex arr_sol[i]
            lnSum = sumAB.complexln();
            arr_sol[i] = new complex(lnSum.b, -lnSum.a); // adding in values to arr_sol
        }
        return arr_sol;
    }
    public static complex[] complexarcsin(complex x){
        return x.complexarcsin();
    }


    /*
     * complexarccos
     * finds the arccosine of a given complex value
     * static and nonstatic methods available
     * static parameter: complex x is the complex that's arccosine will be found
     * returns arr_sol: new complex array object that holds the two solutions of arccossine of a complex value
     */
    public complex[] complexarccos(){ 
        complex xsquare = power(2);
        complex solA = new complex(xsquare.a - 1, xsquare.b);
        complex[] rootA = solA.kroot(2);
        complex[] arr_sol = new complex[2];
        complex sumAB = new complex();
        complex lnSum = new complex();
        for(int i = 0; i < rootA.length; i++){
            sumAB = add(rootA[i]);
            lnSum = sumAB.complexln();
            arr_sol[i] = new complex(lnSum.b, -lnSum.a);
        }
        return arr_sol;
    }
    public static complex[] complexarccos(complex x){
        return x.complexarccos();
    }


    /*
     * complexarctam
     * finds the arctangent of a given complex value
     * static and nonstatic methods available
     * static parameter: complex x is the complex that's arctangent will be found
     * returns sol: new complex object that holds the solution of arctangent of a complex value
     */
    public complex complexarctan(){
        complex solA = new complex(1 + b, -a); // calculations for solA
        solA = solA.complexln();
        complex solB = new complex(1 - b, a); // calculations for solB
        solB = solB.complexln();
        complex solC = solA.subtract(solB); // solC used in calculations for final solution
        complex sol = new complex(0.5 * -solC.b, 0.5 * solC.a); // calculations for sol
        sol.a = (sol.a + Math.PI) % Math.PI; // mutate solA to give a positive real arctangent value
        return sol;
    }
    public static complex complexarctan(complex x){
        return x.complexarctan();
    }


    /*
     * quadratic
     * finds the zeros of a quadratic function given a, b, c are complex values
     * static method available
     * static parameter: complex x is the a value, complex y is the b value, and complex z is the c value, as a,b,c are treated as part of a quadratic ax^2 + bx + c
     * returns arr_sol: new complex array object that holds all zeros of the complex quadratic
     */
    public static complex[] quadratic(complex x, complex y, complex z){
        complex disc1 = y.power(2); // disc1 (short for "discriminant") needed to find the value of the discriminant
        complex disc2 = new complex(4 * x.multiply(z).a, 4 * x.multiply(z).b); // disc2 needed to find the value of the discriminant
        complex disc = disc1.subtract(disc2); // calculations for discriminant
        int discNum = 2; //discriminant number: how many solutions the discriminant yields ( 2 for nonzero values and 1 for zero values)
        if(disc.equals(new complex())){ 
            discNum = 1;
        }
        complex[] arr_sol = new complex[discNum];
        complex negY = new complex(-y.a, -y.b);
        complex den = new complex(2 * x.a, 2 * x.b); // calculations for denominator

        if(disc.equals(new complex())){ // filling in arr_sol if disc equals 0
            complex num = negY.subtract(disc.kroot(2)[0]);
            arr_sol[0] = num.divide(den);
        }
        else{  // filling in arr_sol if disc does not equal 0 (2 solutions)
            complex num = negY.subtract(disc.kroot(2)[0]);
            arr_sol[0] = num.divide(den);
            num = negY.add(disc.kroot(2)[0]);
            arr_sol[1] = num.divide(den);
        }
        return arr_sol;
    }

    
}
