;give wavelength and amplitude for each zeeman component
;Ken Liao
;based on Igor Bespamyatnov's cxrs_fit.pro
;B should be a scalar
;pol is 0 for pi component, -1 or 1 for sigma component.
;But I'm not sure which of -1 or +1 is +circularly polarized, and which is -circularly polarized
pro zeeman_boron,B,wavelength,amplitude,pi_comp_arr,pol
    ;intensity of each line component
    amplitude = [1261828d, 388255d, 1067701d, 16177d, 711801d, 889751d, 48532d, 970637d, 727978d, 97064d, 1164765d, 582382d, 161773d, 1294183d, 452964d, 242659d, 1358892d, 339723d, 339723d, 1358892d, 242659d, 452964d, 1294183d, 161773d,$
    582382d, 1164765d, 97064d, 727978d, 970637d, 48532d, 889751d, 711801d, 16177d, 1067701d, 388255d, 1261828d, 906707d, 329712d, 741851d, 16486d, 593481d, 593481d, 49457d, 791308d, 461596d, 98913d, 923193d, 346197d, 164856d,$
    989135d, 247284d, 247284d, 989135d, 164856d, 346197d, 923193d, 98913d, 461596d, 791308d, 49457d, 593481d, 593481d, 16486d, 741851d, 329712d, 906707d, 32353d, 2941d, 2941d, 21658d, 5348d, 5348d, 13102d, 7219d, 7219d,$
    6685d, 8556d, 8556d, 2406d, 9358d, 9358d, 267d, 9626d, 9626d, 267d, 9358d, 9358d, 2406d, 8556d, 8556d, 6685d, 7219d, 7219d, 13102d, 5348d, 5348d, 21658d, 2941d, 2941d, 32353d, 1242766d, 451915d, 1016809d, 22596d, 813447d,$
    813447d, 67787d, 1084596d, 632681d, 135574d, 1265362d, 474511d, 225957d, 1355745d, 338936d, 338936d, 1355745d, 225957d, 474511d, 1265362d, 135574d, 632681d, 1084596d, 67787d, 813447d, 813447d, 22596d, 1016809d, 451915d,$
    1242766d, 4674d, 1700d, 85d, 3824d, 3059d, 255d, 3059d, 4079d, 510d, 2380d, 4759d, 850d, 1785d, 5099d, 1275d, 1275d, 5099d, 1785d, 850d, 4759d, 2380d, 510d, 4079d, 3059d, 255d, 3059d, 3824d, 85d, 1700d, 4674d, 170d, 19d,$
    19d, 103d, 34d, 34d, 52d, 44d, 44d, 19d, 50d, 50d, 2d, 52d, 52d, 2d, 50d, 50d, 19d, 44d, 44d, 52d, 34d, 34d, 103d, 19d, 19d, 170d, 641163d, 284961d, 498682d, 17810d, 498682d, 374012d, 53430d, 641163d, 267151d, 106861d,$
    712403d, 178101d, 178101d, 712403d, 106861d, 267151d, 641163d, 53430d, 374012d, 498682d, 17810d, 498682d, 284961d, 641163d, 32969d, 3663d, 3663d, 19944d, 6512d, 6512d, 10176d, 8548d, 8548d, 3663d, 9769d, 9769d, 407d,$
    10176d, 10176d, 407d, 9769d, 9769d, 3663d, 8548d, 8548d, 10176d, 6512d, 6512d, 19944d, 3663d, 3663d, 32969d, 886623d, 394054d, 689595d, 24628d, 689595d, 517196d, 73885d, 886623d, 369426d, 147770d, 985136d, 246284d,$
    246284d, 985136d, 147770d, 369426d, 886623d, 73885d, 517196d, 689595d, 24628d, 689595d, 394054d, 886623d, 14708d, 6537d, 409d, 11439d, 11439d, 1226d, 8579d, 14708d, 2451d, 6128d, 16342d, 4085d, 4085d, 16342d, 6128d,$
    2451d, 14708d, 8579d, 1226d, 11439d, 11439d, 409d, 6537d, 14708d, 817d, 117d, 117d, 417d, 200d, 200d, 150d, 250d, 250d, 17d, 267d, 267d, 17d, 250d, 250d, 150d, 200d, 200d, 417d, 117d, 117d, 817d, 441986d, 252563d,$
    315704d, 21047d, 420939d, 210470d, 63141d, 505127d, 126282d, 126282d, 505127d, 63141d, 210470d, 420939d, 21047d, 315704d, 252563d, 441986d, 4570d, 2031d, 127d, 3555d, 3555d, 381d, 2666d, 4570d, 762d, 1904d, 5078d, 1270d,$
    1270d, 5078d, 1904d, 762d, 4570d, 2666d, 381d, 3555d, 3555d, 127d, 2031d, 4570d, 35617d, 5088d, 5088d, 18172d, 8723d, 8723d, 6542d, 10903d, 10903d, 727d, 11630d, 11630d, 727d, 10903d, 10903d, 6542d, 8723d, 8723d, 18172d,$
    5088d, 5088d, 35617d, 618345d, 353340d, 441675d, 29445d, 588900d, 294450d, 88335d, 706680d, 176670d, 176670d, 706680d, 88335d, 294450d, 588900d, 29445d, 441675d, 353340d, 618345d, 30705d, 17545d, 1462d, 21932d, 29242d,$
    4386d, 14621d, 35091d, 8773d, 8773d, 35091d, 14621d, 4386d, 29242d, 21932d, 1462d, 17545d, 30705d, 2925d, 585d, 585d, 1053d, 936d, 936d, 117d, 1053d, 1053d, 117d, 936d, 936d, 1053d, 585d, 585d, 2925d, 290432d, 232346d,$
    174259d, 29043d, 348519d, 87130d, 87130d, 348519d, 29043d, 174259d, 232346d, 290432d, 14182d, 8104d, 675d, 10130d, 13507d, 2026d, 6753d, 16208d, 4052d, 4052d, 16208d, 6753d, 2026d, 13507d, 10130d, 675d, 8104d, 14182d,$
    42088d, 8418d, 8418d, 15152d, 13468d, 13468d, 1684d, 15152d, 15152d, 1684d, 13468d, 13468d, 15152d, 8418d, 8418d, 42088d, 412639d, 330111d, 247583d, 41264d, 495167d, 123792d, 123792d, 495167d, 41264d, 247583d, 330111d,$
    412639d, 52160d, 41728d, 5216d, 31296d, 62593d, 15648d, 15648d, 62593d, 31296d, 5216d, 41728d, 52160d, 10436d, 3479d, 3479d, 1160d, 4638d, 4638d, 1160d, 3479d, 3479d, 10436d, 166995d, 222660d, 55665d, 55665d, 222660d,$
    166995d, 28656d, 22925d, 2866d, 17194d, 34387d, 8597d, 8597d, 34387d, 17194d, 2866d, 22925d, 28656d, 58071d, 19357d, 19357d, 6452d, 25809d, 25809d, 6452d, 19357d, 19357d, 58071d, 242267d, 323023d, 80756d, 80756d, 323023d,$
    242267d, 71798d, 95731d, 23933d, 23933d, 95731d, 71798d, 47926d, 47926d, 47926d, 47926d, 43451d, 57934d, 14484d, 14484d, 57934d, 43451d, 111242d, 111242d, 111242d, 111242d]

    pi_comp_arr = [0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0$
    ,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0$
    ,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0$
    ,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0$
    ,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1$
    ,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1$
    ,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1$
    ,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1$
    ,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1$
    ,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1$
    ,0,1,0,0,1,0,1,0,0,1]

    zeem_coef_1 = [0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, $
    0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d,$
    0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d,$
    0.0002021726252145d, 0.0002021726252145d, 0.0002021726252145d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d,$
    0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d,$
    0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d,$
    0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d,$
    0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d,$
    0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d, 0.0002021700706361d,$
    0.0002021700706361d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d,$
    0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d,$
    0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021756987772d, 0.0002021665182829d, 0.0002021665182829d,$
    0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d,$
    0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d,$
    0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021665182829d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d,$
    0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d,$
    0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d,$
    0.0002021721462261d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d,$
    0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d,$
    0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d,$
    0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d,$
    0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021721462261d, 0.0002021805887286d, 0.0002021805887286d,$
    0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d,$
    0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d, 0.0002021805887286d,$
    0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d,$
    0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d, 0.0002021668375904d,$
    0.0002021668375904d, 0.0002021668375904d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d,$
    0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d, 0.0002021752796495d,$
    0.0002021752796495d, 0.0002021752796495d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d,$
    0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021893313536d, 0.0002021667976769d, 0.0002021667976769d,$
    0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d,$
    0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d, 0.0002021667976769d,$
    0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d,$
    0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d, 0.0002021752397326d,$
    0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d,$
    0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d, 0.0002021893712760d, 0.0002021663985428d, 0.0002021663985428d, 0.0002021663985428d, 0.0002021663985428d,$
    0.0002021663985428d, 0.0002021663985428d, 0.0002021663985428d, 0.0002021663985428d, 0.0002021663985428d, 0.0002021663985428d, 0.0002021663985428d, 0.0002021663985428d, 0.0002021663985428d, 0.0002021663985428d, 0.0002021663985428d,$
    0.0002021663985428d, 0.0002021663985428d, 0.0002021663985428d, 0.0002021805488096d, 0.0002021805488096d, 0.0002021805488096d, 0.0002021805488096d, 0.0002021805488096d, 0.0002021805488096d, 0.0002021805488096d, 0.0002021805488096d,$
    0.0002021805488096d, 0.0002021805488096d, 0.0002021805488096d, 0.0002021805488096d, 0.0002021805488096d, 0.0002021805488096d, 0.0002021805488096d, 0.0002021805488096d, 0.0002022085957413d, 0.0002022085957413d, 0.0002022085957413d,$
    0.0002022085957413d, 0.0002022085957413d, 0.0002022085957413d, 0.0002022085957413d, 0.0002022085957413d, 0.0002022085957413d, 0.0002022085957413d, 0.0002022085957413d, 0.0002022085957413d, 0.0002021663586295d, 0.0002021663586295d,$
    0.0002021663586295d, 0.0002021663586295d, 0.0002021663586295d, 0.0002021663586295d, 0.0002021663586295d, 0.0002021663586295d, 0.0002021663586295d, 0.0002021663586295d, 0.0002021663586295d, 0.0002021663586295d, 0.0002021663586295d,$
    0.0002021663586295d, 0.0002021663586295d, 0.0002021663586295d, 0.0002021663586295d, 0.0002021663586295d, 0.0002021804290529d, 0.0002021804290529d, 0.0002021804290529d, 0.0002021804290529d, 0.0002021804290529d, 0.0002021804290529d,$
    0.0002021804290529d, 0.0002021804290529d, 0.0002021804290529d, 0.0002021804290529d, 0.0002021804290529d, 0.0002021804290529d, 0.0002021804290529d, 0.0002021804290529d, 0.0002021804290529d, 0.0002021804290529d, 0.0002022087155313d,$
    0.0002022087155313d, 0.0002022087155313d, 0.0002022087155313d, 0.0002022087155313d, 0.0002022087155313d, 0.0002022087155313d, 0.0002022087155313d, 0.0002022087155313d, 0.0002022087155313d, 0.0002022087155313d, 0.0002022087155313d,$
    0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d,$
    0.0002021627664948d, 0.0002021910480308d, 0.0002021910480308d, 0.0002021910480308d, 0.0002021910480308d, 0.0002021910480308d, 0.0002021910480308d, 0.0002021910480308d, 0.0002021910480308d, 0.0002021910480308d, 0.0002021910480308d,$
    0.0002022705860938d, 0.0002022705860938d, 0.0002022705860938d, 0.0002022705860938d, 0.0002022705860938d, 0.0002022705860938d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d,$
    0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021627664948d, 0.0002021908284542d, 0.0002021908284542d, 0.0002021908284542d, 0.0002021908284542d,$
    0.0002021908284542d, 0.0002021908284542d, 0.0002021908284542d, 0.0002021908284542d, 0.0002021908284542d, 0.0002021908284542d, 0.0002022755805275d, 0.0002022755805275d, 0.0002022755805275d, 0.0002022755805275d, 0.0002022755805275d,$
    0.0002022755805275d, 0.0002021406976243d, 0.0002021406976243d, 0.0002021406976243d, 0.0002021406976243d, 0.0002021406976243d, 0.0002021406976243d, 0.0002022254675723d, 0.0002022254675723d, 0.0002022254675723d, 0.0002022254675723d,$
    0.0002021376650277d, 0.0002021376650277d, 0.0002021376650277d, 0.0002021376650277d, 0.0002021376650277d, 0.0002021376650277d, 0.0002022171810550d, 0.0002022171810550d, 0.0002022171810550d, 0.0002022171810550d]

    zeem_coef_2 = [-0.0000000046638816d, 0.0000000003587595d, -0.0000000047291106d, 0.0000000053814006d, 0.0000000002935305d, -0.0000000047943396d, 0.0000000053161716d, 0.0000000002283015d, -0.0000000048595686d, 0.0000000052509426d,$
    0.0000000001630725d, -0.0000000049247976d, 0.0000000051857136d, 0.0000000000978435d, -0.0000000049900266d, 0.0000000051204846d, 0.0000000000326145d, -0.0000000050552556d, 0.0000000050552556d, -0.0000000000326145d,$
    -0.0000000051204846d, 0.0000000049900266d, -0.0000000000978435d, -0.0000000051857136d, 0.0000000049247976d, -0.0000000001630725d, -0.0000000052509426d, 0.0000000048595686d, -0.0000000002283015d, -0.0000000053161716d,$
    0.0000000047943396d, -0.0000000002935305d, -0.0000000053814006d, 0.0000000047291106d, -0.0000000003587595d, 0.0000000046638816d, -0.0000000046638805d, 0.0000000004239896d, -0.0000000047581004d, 0.0000000055118597d,$
    0.0000000003297697d, -0.0000000048523204d, 0.0000000054176398d, 0.0000000002355498d, -0.0000000049465403d, 0.0000000053234199d, 0.0000000001413299d, -0.0000000050407602d, 0.0000000052292000d, 0.0000000000471100d,$
    -0.0000000051349801d, 0.0000000051349801d, -0.0000000000471100d, -0.0000000052292000d, 0.0000000050407602d, -0.0000000001413299d, -0.0000000053234199d, 0.0000000049465403d, -0.0000000002355498d, -0.0000000054176398d,$
    0.0000000048523204d, -0.0000000003297697d, -0.0000000055118597d, 0.0000000047581004d, -0.0000000004239896d, 0.0000000046638805d, 0.0000000043051215d, -0.0000000007827487d, 0.0000000086102422d, 0.0000000035223721d,$
    -0.0000000015654980d, 0.0000000078274929d, 0.0000000027396227d, -0.0000000023482474d, 0.0000000070447435d, 0.0000000019568734d, -0.0000000031309967d, 0.0000000062619942d, 0.0000000011741240d, -0.0000000039137461d,$
    0.0000000054792448d, 0.0000000003913747d, -0.0000000046964954d, 0.0000000046964954d, -0.0000000003913747d, -0.0000000054792448d, 0.0000000039137461d, -0.0000000011741240d, -0.0000000062619942d, 0.0000000031309967d,$
    -0.0000000019568734d, -0.0000000070447435d, 0.0000000023482474d, -0.0000000027396227d, -0.0000000078274929d, 0.0000000015654980d, -0.0000000035223721d, -0.0000000086102422d, 0.0000000007827487d, -0.0000000043051215d,$
    -0.0000000045986525d, -0.0000000002935317d, -0.0000000045334232d, 0.0000000040115890d, -0.0000000002283025d, -0.0000000044681940d, 0.0000000040768183d, -0.0000000001630732d, -0.0000000044029647d, 0.0000000041420476d,$
    -0.0000000000978439d, -0.0000000043377354d, 0.0000000042072769d, -0.0000000000326146d, -0.0000000042725061d, 0.0000000042725061d, 0.0000000000326146d, -0.0000000042072769d, 0.0000000043377354d, 0.0000000000978439d,$
    -0.0000000041420476d, 0.0000000044029647d, 0.0000000001630732d, -0.0000000040768183d, 0.0000000044681940d, 0.0000000002283025d, -0.0000000040115890d, 0.0000000045334232d, 0.0000000002935317d, 0.0000000045986525d,$
    0.0000000046638805d, -0.0000000004239896d, -0.0000000055118597d, 0.0000000047581004d, -0.0000000003297697d, -0.0000000054176398d, 0.0000000048523204d, -0.0000000002355498d, -0.0000000053234199d, 0.0000000049465403d,$
    -0.0000000001413299d, -0.0000000052292000d, 0.0000000050407602d, -0.0000000000471100d, -0.0000000051349801d, 0.0000000051349801d, 0.0000000000471100d, -0.0000000050407602d, 0.0000000052292000d, 0.0000000001413299d,$
    -0.0000000049465403d, 0.0000000053234199d, 0.0000000002355498d, -0.0000000048523204d, 0.0000000054176398d, 0.0000000003297697d, -0.0000000047581004d, 0.0000000055118597d, 0.0000000004239896d, -0.0000000046638805d,$
    -0.0000000042398934d, -0.0000000084797849d, 0.0000000009421966d, -0.0000000032976949d, -0.0000000075375864d, 0.0000000018843951d, -0.0000000023554964d, -0.0000000065953878d, 0.0000000028265937d, -0.0000000014132978d,$
    -0.0000000056531893d, 0.0000000037687922d, -0.0000000004710993d, -0.0000000047109908d, 0.0000000047109908d, 0.0000000004710993d, -0.0000000037687922d, 0.0000000056531893d, 0.0000000014132978d, -0.0000000028265937d,$
    0.0000000065953878d, 0.0000000023554964d, -0.0000000018843951d, 0.0000000075375864d, 0.0000000032976949d, -0.0000000009421966d, 0.0000000084797849d, 0.0000000042398934d, -0.0000000046638811d, 0.0000000005182089d,$
    -0.0000000048119408d, 0.0000000057002990d, 0.0000000003701492d, -0.0000000049600005d, 0.0000000055522393d, 0.0000000002220895d, -0.0000000051080602d, 0.0000000054041796d, 0.0000000000740298d, -0.0000000052561199d,$
    0.0000000052561199d, -0.0000000000740298d, -0.0000000054041796d, 0.0000000051080602d, -0.0000000002220895d, -0.0000000055522393d, 0.0000000049600005d, -0.0000000003701492d, -0.0000000057002990d, 0.0000000048119408d,$
    -0.0000000005182089d, 0.0000000046638811d, 0.0000000042398934d, -0.0000000009421966d, 0.0000000084797849d, 0.0000000032976949d, -0.0000000018843951d, 0.0000000075375864d, 0.0000000023554964d, -0.0000000028265937d,$
    0.0000000065953878d, 0.0000000014132978d, -0.0000000037687922d, 0.0000000056531893d, 0.0000000004710993d, -0.0000000047109908d, 0.0000000047109908d, -0.0000000004710993d, -0.0000000056531893d, 0.0000000037687922d,$
    -0.0000000014132978d, -0.0000000065953878d, 0.0000000028265937d, -0.0000000023554964d, -0.0000000075375864d, 0.0000000018843951d, -0.0000000032976949d, -0.0000000084797849d, 0.0000000009421966d, -0.0000000042398934d,$
    -0.0000000045696602d, -0.0000000003297687d, -0.0000000044754406d, 0.0000000039101228d, -0.0000000002355491d, -0.0000000043812209d, 0.0000000040043424d, -0.0000000001413294d, -0.0000000042870013d, 0.0000000040985620d,$
    -0.0000000000471098d, -0.0000000041927817d, 0.0000000041927817d, 0.0000000000471098d, -0.0000000040985620d, 0.0000000042870013d, 0.0000000001413294d, -0.0000000040043424d, 0.0000000043812209d, 0.0000000002355491d,$
    -0.0000000039101228d, 0.0000000044754406d, 0.0000000003297687d, 0.0000000045696602d, 0.0000000046638811d, -0.0000000005182089d, -0.0000000057002990d, 0.0000000048119408d, -0.0000000003701492d, -0.0000000055522393d,$
    0.0000000049600005d, -0.0000000002220895d, -0.0000000054041796d, 0.0000000051080602d, -0.0000000000740298d, -0.0000000052561199d, 0.0000000052561199d, 0.0000000000740298d, -0.0000000051080602d, 0.0000000054041796d,$
    0.0000000002220895d, -0.0000000049600005d, 0.0000000055522393d, 0.0000000003701492d, -0.0000000048119408d, 0.0000000057002990d, 0.0000000005182089d, -0.0000000046638811d, -0.0000000041456726d, -0.0000000082913444d,$
    0.0000000011844772d, -0.0000000029611947d, -0.0000000071068665d, 0.0000000023689551d, -0.0000000017767168d, -0.0000000059223887d, 0.0000000035534329d, -0.0000000005922389d, -0.0000000047379108d, 0.0000000047379108d,$
    0.0000000005922389d, -0.0000000035534329d, 0.0000000059223887d, 0.0000000017767168d, -0.0000000023689551d, 0.0000000071068665d, 0.0000000029611947d, -0.0000000011844772d, 0.0000000082913444d, 0.0000000041456726d,$
    -0.0000000046638811d, 0.0000000006662686d, -0.0000000049303885d, 0.0000000059964184d, 0.0000000003997612d, -0.0000000051968960d, 0.0000000057299109d, 0.0000000001332537d, -0.0000000054634035d, 0.0000000054634035d,$
    -0.0000000001332537d, -0.0000000057299109d, 0.0000000051968960d, -0.0000000003997612d, -0.0000000059964184d, 0.0000000049303885d, -0.0000000006662686d, 0.0000000046638811d, 0.0000000045696602d, 0.0000000003297687d,$
    -0.0000000039101228d, 0.0000000044754406d, 0.0000000002355491d, -0.0000000040043424d, 0.0000000043812209d, 0.0000000001413294d, -0.0000000040985620d, 0.0000000042870013d, 0.0000000000471098d, -0.0000000041927817d,$
    0.0000000041927817d, -0.0000000000471098d, -0.0000000042870013d, 0.0000000040985620d, -0.0000000001413294d, -0.0000000043812209d, 0.0000000040043424d, -0.0000000002355491d, -0.0000000044754406d, 0.0000000039101228d,$
    -0.0000000003297687d, -0.0000000045696602d, 0.0000000041456726d, -0.0000000011844772d, 0.0000000082913444d, 0.0000000029611947d, -0.0000000023689551d, 0.0000000071068665d, 0.0000000017767168d, -0.0000000035534329d,$
    0.0000000059223887d, 0.0000000005922389d, -0.0000000047379108d, 0.0000000047379108d, -0.0000000005922389d, -0.0000000059223887d, 0.0000000035534329d, -0.0000000017767168d, -0.0000000071068665d, 0.0000000023689551d,$
    -0.0000000029611947d, -0.0000000082913444d, 0.0000000011844772d, -0.0000000041456726d, -0.0000000045158211d, -0.0000000003701492d, -0.0000000043677614d, 0.0000000037755226d, -0.0000000002220895d, -0.0000000042197017d,$
    0.0000000039235823d, -0.0000000000740298d, -0.0000000040716420d, 0.0000000040716420d, 0.0000000000740298d, -0.0000000039235823d, 0.0000000042197017d, 0.0000000002220895d, -0.0000000037755226d, 0.0000000043677614d,$
    0.0000000003701492d, 0.0000000045158211d, 0.0000000046638811d, -0.0000000006662686d, -0.0000000059964184d, 0.0000000049303885d, -0.0000000003997612d, -0.0000000057299109d, 0.0000000051968960d, -0.0000000001332537d,$
    -0.0000000054634035d, 0.0000000054634035d, 0.0000000001332537d, -0.0000000051968960d, 0.0000000057299109d, 0.0000000003997612d, -0.0000000049303885d, 0.0000000059964184d, 0.0000000006662686d, -0.0000000046638811d,$
    -0.0000000039976126d, -0.0000000079952247d, 0.0000000015990446d, -0.0000000023985675d, -0.0000000063961797d, 0.0000000031980896d, -0.0000000007995225d, -0.0000000047971347d, 0.0000000047971347d, 0.0000000007995225d,$
    -0.0000000031980896d, 0.0000000063961797d, 0.0000000023985675d, -0.0000000015990446d, 0.0000000079952247d, 0.0000000039976126d, -0.0000000046638811d, 0.0000000009327761d, -0.0000000052857318d, 0.0000000065294333d,$
    0.0000000003109254d, -0.0000000059075826d, 0.0000000059075826d, -0.0000000003109254d, -0.0000000065294333d, 0.0000000052857318d, -0.0000000009327761d, 0.0000000046638811d, 0.0000000045158211d, 0.0000000003701492d,$
    -0.0000000037755226d, 0.0000000043677614d, 0.0000000002220895d, -0.0000000039235823d, 0.0000000042197017d, 0.0000000000740298d, -0.0000000040716420d, 0.0000000040716420d, -0.0000000000740298d, -0.0000000042197017d,$
    0.0000000039235823d, -0.0000000002220895d, -0.0000000043677614d, 0.0000000037755226d, -0.0000000003701492d, -0.0000000045158211d, 0.0000000039976126d, -0.0000000015990446d, 0.0000000079952247d, 0.0000000023985675d,$
    -0.0000000031980896d, 0.0000000063961797d, 0.0000000007995225d, -0.0000000047971347d, 0.0000000047971347d, -0.0000000007995225d, -0.0000000063961797d, 0.0000000031980896d, -0.0000000023985675d, -0.0000000079952247d,$
    0.0000000015990446d, -0.0000000039976126d, -0.0000000043973733d, -0.0000000003997612d, -0.0000000041308659d, 0.0000000035978510d, -0.0000000001332537d, -0.0000000038643584d, 0.0000000038643584d, 0.0000000001332537d,$
    -0.0000000035978510d, 0.0000000041308659d, 0.0000000003997612d, 0.0000000043973733d, 0.0000000046638811d, -0.0000000009327761d, -0.0000000065294333d, 0.0000000052857318d, -0.0000000003109254d, -0.0000000059075826d,$
    0.0000000059075826d, 0.0000000003109254d, -0.0000000052857318d, 0.0000000065294333d, 0.0000000009327761d, -0.0000000046638811d, -0.0000000037311048d, -0.0000000074622095d, 0.0000000024874031d, -0.0000000012437016d,$
    -0.0000000049748063d, 0.0000000049748063d, 0.0000000012437016d, -0.0000000024874031d, 0.0000000074622095d, 0.0000000037311048d, -0.0000000046638811d, 0.0000000015546268d, -0.0000000077731348d, 0.0000000077731348d,$
    -0.0000000015546268d, 0.0000000046638811d, 0.0000000043973733d, 0.0000000003997612d, -0.0000000035978510d, 0.0000000041308659d, 0.0000000001332537d, -0.0000000038643584d, 0.0000000038643584d, -0.0000000001332537d,$
    -0.0000000041308659d, 0.0000000035978510d, -0.0000000003997612d, -0.0000000043973733d, 0.0000000037311048d, -0.0000000024874031d, 0.0000000074622095d, 0.0000000012437016d, -0.0000000049748063d, 0.0000000049748063d,$
    -0.0000000012437016d, -0.0000000074622095d, 0.0000000024874031d, -0.0000000037311048d, -0.0000000040420302d, -0.0000000003109255d, -0.0000000034201792d, 0.0000000034201792d, 0.0000000003109255d, 0.0000000040420302d,$
    0.0000000046638811d, -0.0000000015546268d, -0.0000000077731348d, 0.0000000077731348d, 0.0000000015546268d, -0.0000000046638811d, -0.0000000031092540d, -0.0000000062185077d, 0.0000000062185077d, 0.0000000031092540d,$
    0.0000000040420302d, 0.0000000003109255d, -0.0000000034201792d, 0.0000000034201792d, -0.0000000003109255d, -0.0000000040420302d, 0.0000000031092540d, -0.0000000062185077d, 0.0000000062185077d, -0.0000000031092540d]

    wavelength = 1d0/(zeem_coef_1+B*zeem_coef_2)
    wavelength = wavelength / (1d0 + 2.735182d-4 + 131.4182d/wavelength^2+2.76249d8/wavelength^4) ;Vacuum to air. Morton (1991, ApJS, 77, 119)

    pol = -1*((zeem_coef_2 le 0) and (pi_comp_arr eq 0)) + ((zeem_coef_2 gt 0) and (pi_comp_arr eq 0))
end

pro test_zeeman_boron
    B=5.4
    zeeman_boron,B,wavelength,amplitude,pi_comp_arr,pol
    window,1
    polin=1
    plot,[4941,4947],[0,1.4e6],/nodata
    for i=0,n_elements(wavelength)-1 do begin
	case pol[i] of
	    -1: oplot,[wavelength[i],wavelength[i]],[0,amplitude[i]],color='0000ff'x
	    0: oplot,[wavelength[i],wavelength[i]],[0,amplitude[i]],color='00ff00'x
	    1: oplot,[wavelength[i],wavelength[i]],[0,amplitude[i]],color='ff0000'x
	endcase
    endfor
end
