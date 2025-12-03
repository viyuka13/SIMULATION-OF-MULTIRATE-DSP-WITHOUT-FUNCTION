# EXP 6 : SPEECH RECOGNITION USING SCILAB

## AIM: 

To perform and verify multirate DSP without function using SCILAB.

## APPARATUS REQUIRED: 
PC installed with SCILAB. 

## PROGRAM : 

//  SPEECH RECOGNITION USING SCILAB

```c
// ---------------------------
// Multirate DSP demonstration (no user-defined functions)
// - generate signal
// - upsample by L (zero insertion)
// - low-pass filter (sinc FIR with Hamming window)
// - decimate by M
// - show time & frequency-domain results
// ---------------------------

// Parameters
fs = 1000;                // original sampling freq (Hz)
T  = 0.1;                 // signal duration (s)
t  = 0:1/fs:T-1/fs;       // time vector

// Test signal: two sinusoids (50 Hz and 120 Hz) + small noise
x = sin(2*%pi*50*t) + 0.5*sin(2*%pi*120*t) + 0.05*rand(1,length(t));

// Multirate factors
L = 5;   // up-sampling factor (interpolation)
M = 8;   // down-sampling factor (decimation)

// --- 1) Upsample by L (zero insertion) ---
Nx = length(x);
xu = zeros(1, L*Nx);                    // zero-stuffed vector
xu(1:L:L*Nx) = x;                       // place original samples every L-th position

// New sampling rate after upsampling
fs_up = fs * L;

// --- 2) Design lowpass FIR filter (sinc windowed by Hamming)
//     We design it directly (no external function calls).
N = 101;                                 // filter length (odd preferred)
n = 0:N-1;
center = (N-1)/2;

// cutoff frequency (Hz) for interpolation filter: keep original Nyquist (fs/2)
// normalized to new sampling rate fs_up:
fcut = (fs/2);                           // desired cutoff in Hz (preserve original Nyquist)
fcut_norm = fcut / fs_up;                // normalized (cycles/sample at fs_up) -> between 0..0.5

// Create sinc kernel: h[n] = 2*fcut_norm * sinc(2*fcut_norm*(n-center))
// sinc(z) = sin(pi*z) / (pi*z)
xarg = 2*fcut_norm * (n - center);
h = zeros(1,N);
epsv = 1e-12;
for k=1:N
    if abs(xarg(k)) < epsv then
        s = 1;
    else
        s = sin(%pi * xarg(k)) / (%pi * xarg(k));
    end
    h(k) = 2 * fcut_norm * s;
end

// Apply Hamming window
w = 0.54 - 0.46 * cos(2*%pi*n/(N-1));
h = h .* w;

// Normalize filter gain at DC (optional)
h = h / sum(h);

// --- 3) Filter the upsampled signal (convolution)
y_up = conv(xu, h);            // length = length(xu)+N-1

// Compensate filter delay (group delay = (N-1)/2 samples at fs_up)
delay = floor((N-1)/2);
y_up = y_up(delay+1 : delay+Nx*L);  // remove transient so length matches L*Nx (keeps core part)

// --- 4) Now downsample by M (decimation)
y_dec = y_up(1:M:length(y_up));     // pick every M-th sample
fs_out = fs_up / M;                 // output sampling rate

// --- 5) Plot time-domain signals
scf(0);
clf();
subplot(3,1,1);
plot(t, x);
xtitle("Original signal (time domain)");
xlabel("Time (s)");
ylabel("Amplitude");

subplot(3,1,2);
tu = 0:1/fs_up:(length(xu)-1)/fs_up;
plot(tu, xu);
xtitle("Upsampled (zero-stuffed) signal (time domain)");
xlabel("Time (s)");
ylabel("Amplitude");

subplot(3,1,3);
tout = 0:1/fs_out:(length(y_dec)-1)/fs_out;
plot(tout, y_dec);
xtitle("After interpolation (filter) and decimation (time domain)");
xlabel("Time (s)");
ylabel("Amplitude");

// --- 6) Frequency domain: compute magnitude spectra
nfft = 4096;
X  = abs(fft(x, nfft));
X  = X(1:nfft/2) / max(X);
Xu = abs(fft(xu, nfft));
Xu = Xu(1:nfft/2) / max(Xu);
Yd = abs(fft(y_dec, nfft));
Yd = Yd(1:nfft/2) / max(Yd);

f_x  = (0:nfft/2-1) * (fs/nfft);
f_xu = (0:nfft/2-1) * (fs_up/nfft);
f_yd = (0:nfft/2-1) * (fs_out/nfft);

scf(1);
clf();
subplot(3,1,1);
plot(f_x, 20*log10(X+1e-12));
xgrid();
xtitle("Original signal spectrum");
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
xaxis([0 fs/2]);

subplot(3,1,2);
plot(f_xu, 20*log10(Xu+1e-12));
xgrid();
xtitle("Upsampled spectrum (shows images at multiples of fs)");
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
xaxis([0 fs_up/2]);

subplot(3,1,3);
plot(f_yd, 20*log10(Yd+1e-12));
xgrid();
xtitle("Output (after filtering + decimation) spectrum");
xlabel("Frequency (Hz)");
ylabel("Magnitude (dB)");
xaxis([0 fs_out/2]);

// --- 7) Print some diagnostics
mprintf("\nOriginal fs = %g Hz  | Upsampled fs = %g Hz | Output fs = %g Hz\n", fs, fs_up, fs_out);
mprintf("Original length = %d samples | After upsample = %d | After decimate = %d\n", Nx, length(xu), length(y_dec));
mprintf("Filter length N = %d | cutoff = %g Hz (normalized = %g)\n", N, fcut, fcut_norm);
```

## OUTPUT: 

<img width="795" height="419" alt="6" src="https://github.com/user-attachments/assets/ad64b23e-7b5e-4542-80b6-08e8ce661f71" />


## RESULT: 

Thus the decimation process by a factor M and interpolation process by a factor L using 
SCILAB was implemented. 
