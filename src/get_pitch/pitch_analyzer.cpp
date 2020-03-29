/// @file

#include <iostream>
#include <math.h>
#include "pitch_analyzer.h"

using namespace std;

/// Name space of UPC
namespace upc {
  void PitchAnalyzer::autocorrelation(const vector<float> &x, vector<float> &r) const {

    for (unsigned int l = 0; l < r.size(); ++l) {
  		/// \HECHO Compute the autocorrelation r[l]
      for (unsigned int k = l; k < r.size(); ++k){
        r[l] = 0;
        r[l] = r[l] + x[k-l]*x[k];
      }
      r[l] /= r.size();
    }
    
    if (r[0] == 0.0F) //to avoid log() and divide zero 
      r[0] = 1e-10; 
  }

  void PitchAnalyzer::set_window(Window win_type) {
    if (frameLen == 0)
      return;

    window.resize(frameLen);

    switch (win_type) {
    case HAMMING:
      /// \HECHO Implement the Hamming window
      for (unsigned int i = 0; i < frameLen; ++i){
        window[i] = 25/46 - (1-25/46)*cos((2*M_PI*i)/(frameLen-1)); 
      }

      break;
    case RECT:
    default:
      window.assign(frameLen, 1);
    }
  }

  void PitchAnalyzer::set_f0_range(float min_F0, float max_F0) {
    npitch_min = (unsigned int) samplingFreq/max_F0;
    if (npitch_min < 2)
      npitch_min = 2;  // samplingFreq/2

    npitch_max = 1 + (unsigned int) samplingFreq/min_F0;

    //frameLen should include at least 2*T0
    if (npitch_max > frameLen/2)
      npitch_max = frameLen/2;
  }

  bool PitchAnalyzer::unvoiced(float pot, float r1norm, float rmaxnorm, vector<float> x) const {
    /// \TODO Implement a rule to decide whether the sound is voiced or not.
    /// * You can use the standard features (pot, r1norm, rmaxnorm),
    ///   or compute and use other ones.
    float zcr = compute_zeros(x);
    if((r1norm < 0.8 || rmaxnorm < 0.3) && zcr < 400){ //aparentemente la tasa de cruces por 0 aumenta en plan un montón como de 100 a 800, pos más o menos 400, pero busquemos si es posible un valor experimental
      return true;
    } else {
      return false;
    }
    /// \HECHO Agregar el uso de tasa de cruces por 0 a la detección de unvoiced
  }

  float PitchAnalyzer::compute_pitch(vector<float> & x) const {
    if (x.size() != frameLen)
      return -1.0F;

    //Window input frame
    for (unsigned int i=0; i<x.size(); ++i)
      x[i] *= window[i];

    vector<float> r(npitch_max);

    //Compute correlation
    autocorrelation(x, r);

    vector<float>::const_iterator iR = r.begin(), iRMax = iR, iRMin = iR;

    /// \TODO 
	/// Find the lag of the maximum value of the autocorrelation away from the origin.<br>
	/// Choices to set the minimum value of the lag are:
	///    - The first negative value of the autocorrelation.
	///    - The lag corresponding to the maximum value of the pitch.
    ///	   .
	/// In either case, the lag should not exceed that of the minimum value of the pitch.

    for(iR+=npitch_min; iR < r.begin()+npitch_max; ++iR){
      if(*iR > *iRMax){
        iRMax = iR;
      }
      if(r.begin() == iRMin && *iR < 0){
        iRMin = iR;
      }
    }
    
    unsigned int lag = iRMax - r.begin();

    /*if(lag == 0){
      lag= iRMin - r.begin();
    }*/
    
    float pot = 10 * log10(r[0]);

    //You can print these (and other) features, look at them using wavesurfer
    //Based on that, implement a rule for unvoiced
    //change to #if 1 and compile
#if 0
    if (r[0] > 0.0F)
      cout << pot << '\t' << r[1]/r[0] << '\t' << r[lag]/r[0] << endl;
#endif
    
    if (unvoiced(pot, r[1]/r[0], r[lag]/r[0], x))
      return 0;
    else
      return (float) samplingFreq/(float) lag;
  }

  float PitchAnalyzer::compute_zeros(vector<float> x) const {
    float zeros = 0;
    for(int i = 0; i < x.size() - 1; i++) {
      if (x[i] * x[i + 1] > 0) {
        zeros++;
      }
    } 
    return samplingFreq*zeros/(2*(x.size() - 1));
    /// \HECHO Función de tasa de cruces por 0 para detectar voiced sound
  }
  int PitchAnalyzer::compute_lag(vector<float> r) const {
    int aux, lag;
    float max = 0;
    for(int i = 0; i < r.size(); i++) {
      if (r[i] < 0) {
        aux = i;
      }
    }
    // Buscamos el primer valor negativo tal que evitamos el pico del origen
    for(int i = 0; i < r.size() - aux; i++) {
      if (r[i] > max) {
        max = r[i];
        lag = i;
      }
    }
    return lag;
    /// \HECHO Función de encontrar la posición del siguiente pico fuera del origen
  }

  void PitchAnalyzer::median_filter(vector<float> &r) {
    vector<float> auxVector = r;
    vector<float> window(3);
    int j = 0;
    for(int i = 0; i < r.size() - 2; i++) {
      window[0] = auxVector[i], window[1] = auxVector[i+1], window[2] = auxVector[i+2];
      while(!(window[1] > window[0] && window[1] < window[2])) {
        if (window[j%3] > window[(j+1)%3]) {
          swap(window[j%3], window[(j+1)%3]);
        }
        j++;
      }
      r[i] = window[1];
      j = 0;
    }
  }

  void PitchAnalyzer::swap(float &a, float &b) {
    a = a + b;
    b = a - b;
    a = a - b;
  }
}
