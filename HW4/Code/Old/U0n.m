function result = U0n(x)

%     if x < 0.3
%         error('ERROR.\nInput to U0n must be greater than 0.3.\nInput value was %s.', round(x, 4));
%     end

    result = 0.1493 + 0.8493*exp(-61.79*x) + 0.3824*exp(-665.8*x) ...
            - exp(39.42*x - 41.92) - 0.03131*atan(25.59*x - 4.099) ...
            - 0.009434*atan(32.49*x - 15.74);

end