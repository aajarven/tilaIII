float harmonic_kahan(int N){
    float s = 0;
    float x, y, t;
    float e = 0;

    for (int i=1; i<=N; i++){
        x = 1.0/i;
        y = x-e;
        t = s+y;
        e = (t-s)-y;
        s = t;
    }
    
    return s;
}
