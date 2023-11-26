index = -1; best = inf;

for(ii = 0; ii < L; ++ii) {
   vc[ii] = !vc[ii];
   val = eval(vc);
   
   if(val < best) {
     best = val;
     index = ii;
   }
   vc[ii] = !vc[ii];
}

if(best < eval(vc)) vc[index] = !vc[index];