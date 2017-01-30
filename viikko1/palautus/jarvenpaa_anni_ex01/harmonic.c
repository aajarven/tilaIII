float harmonic(void){
    float sum = 0;
    float previousSum;
    int k = 0;

    do {
        k++;
        previousSum = sum;
        sum += 1.0/k;
    } while (sum != previousSum);

    return sum;
}


float harmonic_bunch(int N){
    float sum = 0;
    float previousSum;
    int k = 0;

    do {
        previousSum = sum;
        float bunchSum = 0;

        for (int i=1; i<=N; i++){
            bunchSum += 1.0/(k+i);
        }
        k += N;

        sum += bunchSum;
    } while (sum != previousSum);

    return sum;
}
