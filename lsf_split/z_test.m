% dup = find(diff(index)==0);
x = wavread('../data/source_down/t01s004540.wav');       % source

e_x = stenergy(x);

s_x = strip_sil(x);

figure
subplot(311);
plot(x)
subplot(312);
plot(e_x);
subplot(313);
plot(s_x);