
%Estimate wavelet coeffs at each scale
%pass in scale and compute coherence

function wy=process_wave(WAVE_scale)

%[WAVE, scale]=wavelet_comp(X);

%for scale_no=1:length(scale)
pts_rem=12;
LEN=length(WAVE_scale);
n=LEN-2*pts_rem;
new_WAVE=WAVE_scale(pts_rem:LEN-pts_rem);

for i=1:n
    real_wave=real(new_WAVE(i));
    imag_wave=imag(new_WAVE(i));
    wy(i)=new_WAVE(i);
    
end

end
