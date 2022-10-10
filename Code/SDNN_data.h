// File including the neural network data


#ifndef SDNN_DATA_H
#define SDNN_DATA_H


#include <array>


constexpr std::array<double, 112> W1_data =
{
    -0.935290455818176,
    -0.500665247440338,
    0.655252516269684,
    0.539143681526184,
    0.326104819774628,
    0.0197184197604656,
    -0.424599796533585,
    1.25528955459595,
    1.06807947158813,
    0.177635446190834,
    -0.94041383266449,
    -0.621369600296021,
    0.489421755075455,
    -0.880155622959137,
    -0.693257689476013,
    -0.768862664699554,
    -0.407252311706543,
    0.0863954797387123,
    -4.19955825805664,
    -0.788989901542664,
    2.16184043884277,
    0.950272798538208,
    0.955055236816406,
    1.35560476779938,
    -1.53553128242493,
    -0.722421944141388,
    1.016352891922,
    2.05182528495789,
    -1.30900967121124,
    2.09610962867737,
    0.923514902591705,
    -0.426428288221359,
    -0.444702118635178,
    -1.97467017173767,
    1.77588725090027,
    0.775159418582916,
    0.453386962413788,
    0.979997038841248,
    -1.69584560394287,
    0.629230797290802,
    0.431812435388565,
    -1.74968636035919,
    -0.819498300552368,
    -0.377476781606674,
    -0.0807713270187378,
    -2.4555823802948,
    0.411506712436676,
    0.076002225279808,
    -1.4065614938736,
    -1.32676947116852,
    1.74782061576843,
    0.663716554641724,
    -0.666523456573486,
    -1.35695481300354,
    2.31532955169678,
    1.82193696498871,
    2.23317861557007,
    2.10953140258789,
    -1.26467871665955,
    1.78028631210327,
    -1.88005316257477,
    -0.784120857715607,
    -1.84410011768341,
    -1.0054612159729,
    -0.541405916213989,
    2.65736794471741,
    0.242334827780724,
    0.275797545909882,
    1.33965587615967,
    -0.0816840082406998,
    -3.1013650894165,
    1.88293039798737,
    0.665224373340607,
    -1.4354259967804,
    -0.119827702641487,
    0.245910555124283,
    -0.37157815694809,
    -0.552303791046143,
    0.617477238178253,
    -0.32582214474678,
    0.468512326478958,
    1.16321992874146,
    3.04002141952515,
    0.0838300734758377,
    -0.570439100265503,
    2.18224024772644,
    1.50094819068909,
    -1.69573533535004,
    -0.93999171257019,
    -0.606016755104065,
    -2.35580396652222,
    0.857528805732727,
    -0.908167123794556,
    0.505348861217499,
    0.310234934091568,
    1.70248878002167,
    -0.323960304260254,
    -0.631981253623962,
    0.745268642902374,
    0.200184732675552,
    0.370018124580383,
    0.103307791054249,
    -0.472307235002518,
    1.02071762084961,
    1.36696887016296,
    0.173805519938469,
    -0.18342013657093,
    -0.972113132476807,
    0.465909570455551,
    -0.244284465909004,
    -0.487185537815094,
    -0.514231860637665   



// -0.935290455818176,
// -0.407252311706543,
// -0.444702118635178,
// -1.4065614938736,
// -0.541405916213989,
// 0.468512326478958,
// -0.323960304260254,
// -0.500665247440338,
// 0.0863954797387123,
// -1.97467017173767,
// -1.32676947116852,
// 2.65736794471741,
// 1.16321992874146,
// -0.631981253623962,
// 0.655252516269684,
// -4.19955825805664,
// 1.77588725090027,
// 1.74782061576843,
// 0.242334827780724,
// 3.04002141952515,
// 0.745268642902374,
// 0.539143681526184,
// -0.788989901542664,
// 0.775159418582916,
// 0.663716554641724,
// 0.275797545909882,
// 0.0838300734758377,
// 0.200184732675552,
// 0.326104819774628,
// 2.16184043884277,
// 0.453386962413788,
// -0.666523456573486,
// 1.33965587615967,
// -0.570439100265503,
// 0.370018124580383,
// 0.0197184197604656,
// 0.950272798538208,
// 0.979997038841248,
// -1.35695481300354,
// -0.0816840082406998,
// 2.18224024772644,
// 0.103307791054249,
// -0.424599796533585,
// 0.955055236816406,
// -1.69584560394287,
// 2.31532955169678,
// -3.1013650894165,
// 1.50094819068909,
// -0.472307235002518,
// 1.25528955459595,
// 1.35560476779938,
// 0.629230797290802,
// 1.82193696498871,
// 1.88293039798737,
// -1.69573533535004,
// 1.02071762084961,
// 1.06807947158813,
// -1.53553128242493,
// 0.431812435388565,
// 2.23317861557007,
// 0.665224373340607,
// -0.93999171257019,
// 1.36696887016296,
// 0.177635446190834,
// -0.722421944141388,
// -1.74968636035919,
// 2.10953140258789,
// -1.4354259967804,
// -0.606016755104065,
// 0.173805519938469,
// -0.94041383266449,
// 1.016352891922,
// -0.819498300552368,
// -1.26467871665955,
// -0.119827702641487,
// -2.35580396652222,
// -0.18342013657093,
// -0.621369600296021,
// 2.05182528495789,
// -0.377476781606674,
// 1.78028631210327,
// 0.245910555124283,
// 0.857528805732727,
// -0.972113132476807,
// 0.489421755075455,
// -1.30900967121124,
// -0.0807713270187378,
// -1.88005316257477,
// -0.37157815694809,
// -0.908167123794556,
// 0.465909570455551,
// -0.880155622959137,
// 2.09610962867737,
// -2.4555823802948,
// -0.784120857715607,
// -0.552303791046143,
// 0.505348861217499,
// -0.244284465909004,
// -0.693257689476013,
// 0.923514902591705,
// 0.411506712436676,
// -1.84410011768341,
// 0.617477238178253,
// 0.310234934091568,
// -0.487185537815094,
// -0.768862664699554,
// -0.426428288221359,
// 0.076002225279808,
// -1.0054612159729,
// -0.32582214474678,
// 1.70248878002167,
// -0.514231860637665 
};

constexpr std::array<double, 256> W2_data =
{
0.360904276371002,
0.799066245555878,
-0.137890473008156,
-0.553534328937531,
-0.26813217997551,
0.410066068172455,
0.673012435436249,
-0.0291767455637455,
0.318111479282379,
0.0533493347465992,
-0.508570730686188,
-1.25049459934235,
-0.350822597742081,
0.503553807735443,
0.338341265916824,
-0.382557362318039,
0.588403880596161,
0.892927348613739,
0.45459058880806,
0.151375442743301,
-0.114957578480244,
2.05626368522644,
0.0304539557546377,
-0.365889042615891,
-1.23467755317688,
-1.24699985980988,
-1.39237058162689,
-0.276244521141052,
-0.0544555075466633,
-0.088409535586834,
0.104060463607311,
-1.70834279060364,
-0.238414019346237,
-0.0805100798606873,
0.963208854198456,
-0.213436275720596,
0.495501607656479,
-1.31731188297272,
0.598147571086884,
-0.547534227371216,
1.92111241817474,
-0.0799884423613548,
0.233881056308746,
-1.04845428466797,
-0.426365792751312,
-1.08593332767487,
0.254531353712082,
0.243203043937683,
0.713051199913025,
-0.185351267457008,
0.0768399685621262,
0.374074876308441,
0.0290460418909788,
-0.3818499147892,
0.0576348900794983,
0.44148576259613,
0.239213019609451,
0.249190703034401,
0.222726196050644,
0.767264127731323,
0.164795368909836,
-0.00492732832208276,
-0.399862855672836,
0.206261858344078,
-0.649576544761658,
-0.46263974905014,
-0.841208577156067,
0.40293762087822,
0.0983191505074501,
1.13242828845978,
-0.548608958721161,
0.78550112247467,
0.199957802891731,
0.012489065527916,
-1.15884244441986,
-0.851579666137695,
0.0239765383303165,
0.122826352715492,
-0.393006831407547,
1.11764514446259,
-0.325343698263168,
0.0633938014507294,
-0.454031765460968,
0.293165773153305,
-0.452307194471359,
-0.194093555212021,
-0.828351616859436,
-0.274975776672363,
0.47686380147934,
-0.628260016441345,
-0.0327724367380142,
0.202997699379921,
0.316445022821426,
1.60298192501068,
0.354006260633469,
-0.646530270576477,
0.867284893989563,
-0.145140558481216,
1.720210313797,
-1.2731865644455,
-0.855415821075439,
-0.286877512931824,
-0.707320332527161,
-1.0063316822052,
-1.46386623382568,
1.8688633441925,
0.706775009632111,
0.564194738864899,
-0.257487535476685,
-0.337442010641098,
2.02246809005737,
-0.3382867872715,
0.482344567775726,
0.477101355791092,
-1.09962618350983,
-0.0915534645318985,
-0.913465321063995,
1.11050164699554,
0.0647257938981056,
0.591895937919617,
0.82884556055069,
0.154344439506531,
0.403138846158981,
-0.434597492218018,
0.0867561921477318,
-1.81354510784149,
0.838797748088837,
1.91315817832947,
0.7725670337677,
0.182450532913208,
-0.209928885102272,
0.833817958831787,
0.0271809305995703,
-0.155889749526978,
-0.0105054648593068,
0.824482917785645,
-0.075316034257412,
-0.643773794174194,
1.23831009864807,
2.02473068237305,
-0.185569003224373,
-1.55027294158936,
-0.333671897649765,
-0.112014256417751,
0.360565662384033,
0.879886746406555,
1.37325954437256,
-0.0635737255215645,
-0.566502809524536,
-0.316234946250916,
-0.472098916769028,
0.287602663040161,
-0.965478003025055,
0.289166867733002,
1.45490074157715,
1.17455720901489,
0.0509644858539104,
-1.38335919380188,
1.08752262592316,
-0.895048439502716,
-1.59024369716644,
0.693190693855286,
-0.184622436761856,
-0.39596688747406,
0.0819838121533394,
0.0328229106962681,
1.02089607715607,
-0.443067908287048,
-0.0521412305533886,
0.154734209179878,
-0.177607133984566,
-0.627282559871674,
-0.654224932193756,
-0.416472226381302,
0.0527433417737484,
0.522484421730042,
-0.742458820343018,
-0.736514806747437,
0.512675106525421,
-0.43732076883316,
1.06743907928467,
-0.410609751939774,
-1.63111174106598,
0.482816517353058,
0.0158156920224428,
0.513688683509827,
-0.649710416793823,
-1.2585860490799,
0.815254807472229,
-0.146136447787285,
-1.24378550052643,
0.593237996101379,
0.545360863208771,
-0.149563565850258,
1.43794238567352,
1.25237727165222,
1.66475164890289,
-0.836300373077393,
0.263539284467697,
0.794520556926727,
-1.13636207580566,
0.539906919002533,
0.516298949718475,
-1.34203398227692,
-0.0991810336709023,
-0.746416628360748,
-1.09837651252747,
0.0639115199446678,
-1.58861362934113,
0.101251758635044,
0.942973554134369,
-0.371088027954102,
0.0304836891591549,
-0.11656192690134,
-0.231096774339676,
-0.210966527462006,
-0.792537689208984,
-0.314360707998276,
-0.910883724689484,
-0.792623996734619,
0.261203020811081,
0.20776230096817,
1.16143023967743,
-0.5188929438591,
-0.433489352464676,
-0.462003767490387,
-0.28238382935524,
-0.0114116631448269,
0.462984383106232,
0.299591153860092,
0.412672191858292,
-0.212000951170921,
-0.837797820568085,
0.715392351150513,
-1.66701686382294,
-0.659776329994202,
0.150112211704254,
1.88140320777893,
-0.21598893404007,
0.588217794895172,
1.02075850963593,
0.0726983472704887,
0.628080606460571,
-0.252160668373108,
0.0651043653488159,
0.0916410386562347,
-0.151849120855331,
-0.225768059492111,
0.225631862878799,
0.729324340820313,
-0.389605432748795,
-0.699518501758575,
0.393852144479752,
1.17116713523865,
0.221739679574966,
-1.16239178180695



// 0.360904276371002,
// 0.588403880596161,
// -0.238414019346237,
// 0.713051199913025,
// -0.649576544761658,
// -0.325343698263168,
// 0.867284893989563,
// 0.482344567775726,
// 0.7725670337677,
// 0.360565662384033,
// -1.59024369716644,
// -0.742458820343018,
// 0.545360863208771,
// -1.58861362934113,
// -0.433489352464676,
// 1.02075850963593,
// 0.799066245555878,
// 0.892927348613739,
// -0.0805100798606873,
// -0.185351267457008,
// -0.46263974905014,
// 0.0633938014507294,
// -0.145140558481216,
// 0.477101355791092,
// 0.182450532913208,
// 0.879886746406555,
// 0.693190693855286,
// -0.736514806747437,
// -0.149563565850258,
// 0.101251758635044,
// -0.462003767490387,
// 0.0726983472704887,
// -0.137890473008156,
// 0.45459058880806,
// 0.963208854198456,
// 0.0768399685621262,
// -0.841208577156067,
// -0.454031765460968,
// 1.720210313797,
// -1.09962618350983,
// -0.209928885102272,
// 1.37325954437256,
// -0.184622436761856,
// 0.512675106525421,
// 1.43794238567352,
// 0.942973554134369,
// -0.28238382935524,
// 0.628080606460571,
// -0.553534328937531,
// 0.151375442743301,
// -0.213436275720596,
// 0.374074876308441,
// 0.40293762087822,
// 0.293165773153305,
// -1.2731865644455,
// -0.0915534645318985,
// 0.833817958831787,
// -0.0635737255215645,
// -0.39596688747406,
// -0.43732076883316,
// 1.25237727165222,
// -0.371088027954102,
// -0.0114116631448269,
// -0.252160668373108,
// -0.26813217997551,
// -0.114957578480244,
// 0.495501607656479,
// 0.0290460418909788,
// 0.0983191505074501,
// -0.452307194471359,
// -0.855415821075439,
// -0.913465321063995,
// 0.0271809305995703,
// -0.566502809524536,
// 0.0819838121533394,
// 1.06743907928467,
// 1.66475164890289,
// 0.0304836891591549,
// 0.462984383106232,
// 0.0651043653488159,
// 0.410066068172455,
// 2.05626368522644,
// -1.31731188297272,
// -0.3818499147892,
// 1.13242828845978,
// -0.194093555212021,
// -0.286877512931824,
// 1.11050164699554,
// -0.155889749526978,
// -0.316234946250916,
// 0.0328229106962681,
// -0.410609751939774,
// -0.836300373077393,
// -0.11656192690134,
// 0.299591153860092,
// 0.0916410386562347,
// 0.673012435436249,
// 0.0304539557546377,
// 0.598147571086884,
// 0.0576348900794983,
// -0.548608958721161,
// -0.828351616859436,
// -0.707320332527161,
// 0.0647257938981056,
// -0.0105054648593068,
// -0.472098916769028,
// 1.02089607715607,
// -1.63111174106598,
// 0.263539284467697,
// -0.231096774339676,
// 0.412672191858292,
// -0.151849120855331,
// -0.0291767455637455,
// -0.365889042615891,
// -0.547534227371216,
// 0.44148576259613,
// 0.78550112247467,
// -0.274975776672363,
// -1.0063316822052,
// 0.591895937919617,
// 0.824482917785645,
// 0.287602663040161,
// -0.443067908287048,
// 0.482816517353058,
// 0.794520556926727,
// -0.210966527462006,
// -0.212000951170921,
// -0.225768059492111,
// 0.318111479282379,
// -1.23467755317688,
// 1.92111241817474,
// 0.239213019609451,
// 0.199957802891731,
// 0.47686380147934,
// -1.46386623382568,
// 0.82884556055069,
// -0.075316034257412,
// -0.965478003025055,
// -0.0521412305533886,
// 0.0158156920224428,
// -1.13636207580566,
// -0.792537689208984,
// -0.837797820568085,
// 0.225631862878799,
// 0.0533493347465992,
// -1.24699985980988,
// -0.0799884423613548,
// 0.249190703034401,
// 0.012489065527916,
// -0.628260016441345,
// 1.8688633441925,
// 0.154344439506531,
// -0.643773794174194,
// 0.289166867733002,
// 0.154734209179878,
// 0.513688683509827,
// 0.539906919002533,
// -0.314360707998276,
// 0.715392351150513,
// 0.729324340820313,
// -0.508570730686188,
// -1.39237058162689,
// 0.233881056308746,
// 0.222726196050644,
// -1.15884244441986,
// -0.0327724367380142,
// 0.706775009632111,
// 0.403138846158981,
// 1.23831009864807,
// 1.45490074157715,
// -0.177607133984566,
// -0.649710416793823,
// 0.516298949718475,
// -0.910883724689484,
// -1.66701686382294,
// -0.389605432748795,
// -1.25049459934235,
// -0.276244521141052,
// -1.04845428466797,
// 0.767264127731323,
// -0.851579666137695,
// 0.202997699379921,
// 0.564194738864899,
// -0.434597492218018,
// 2.02473068237305,
// 1.17455720901489,
// -0.627282559871674,
// -1.2585860490799,
// -1.34203398227692,
// -0.792623996734619,
// -0.659776329994202,
// -0.699518501758575,
// -0.350822597742081,
// -0.0544555075466633,
// -0.426365792751312,
// 0.164795368909836,
// 0.0239765383303165,
// 0.316445022821426,
// -0.257487535476685,
// 0.0867561921477318,
// -0.185569003224373,
// 0.0509644858539104,
// -0.654224932193756,
// 0.815254807472229,
// -0.0991810336709023,
// 0.261203020811081,
// 0.150112211704254,
// 0.393852144479752,
// 0.503553807735443,
// -0.088409535586834,
// -1.08593332767487,
// -0.00492732832208276,
// 0.122826352715492,
// 1.60298192501068,
// -0.337442010641098,
// -1.81354510784149,
// -1.55027294158936,
// -1.38335919380188,
// -0.416472226381302,
// -0.146136447787285,
// -0.746416628360748,
// 0.20776230096817,
// 1.88140320777893,
// 1.17116713523865,
// 0.338341265916824,
// 0.104060463607311,
// 0.254531353712082,
// -0.399862855672836,
// -0.393006831407547,
// 0.354006260633469,
// 2.02246809005737,
// 0.838797748088837,
// -0.333671897649765,
// 1.08752262592316,
// 0.0527433417737484,
// -1.24378550052643,
// -1.09837651252747,
// 1.16143023967743,
// -0.21598893404007,
// 0.221739679574966,
// -0.382557362318039,
// -1.70834279060364,
// 0.243203043937683,
// 0.206261858344078,
// 1.11764514446259,
// -0.646530270576477,
// -0.3382867872715,
// 1.91315817832947,
// -0.112014256417751,
// -0.895048439502716,
// 0.522484421730042,
// 0.593237996101379,
// 0.0639115199446678,
// -0.5188929438591,
// 0.588217794895172,
// -1.16239178180695
};


constexpr std::array<double, 256> W3_data =
{
    0.928747475147247,
    0.118580423295498,
    0.0733854696154594,
    -0.216320112347603,
    0.169354125857353,
    -0.428392708301544,
    0.843611776828766,
    0.00536061264574528,
    -0.137300744652748,
    -0.49064639210701,
    0.230988517403603,
    -1.13217556476593,
    0.981262445449829,
    -0.275694042444229,
    0.458381682634354,
    -0.859836757183075,
    -0.22091244161129,
    -0.443158328533173,
    -0.997133553028107,
    -0.994887948036194,
    -0.131350949406624,
    -0.341392010450363,
    0.456412494182587,
    0.00402470072731376,
    0.350713819265366,
    0.262174665927887,
    -0.048433355987072,
    0.543828129768372,
    0.915902376174927,
    0.681087553501129,
    -0.425360172986984,
    0.208959579467773,
    0.340678840875626,
    1.26924169063568,
    0.763332605361938,
    1.23038351535797,
    0.810515403747559,
    0.351482480764389,
    0.302027016878128,
    -0.450287461280823,
    1.04432547092438,
    0.673821926116943,
    -0.16774320602417,
    1.87203669548035,
    0.0531666539609432,
    0.374954998493195,
    0.566317498683929,
    0.792478382587433,
    0.232478484511375,
    0.27703258395195,
    0.256476700305939,
    -0.456453740596771,
    0.390097618103027,
    -0.0704145431518555,
    -1.4993953704834,
    -0.247053876519203,
    -0.0885710567235947,
    1.01270878314972,
    -0.245661944150925,
    -0.478165984153748,
    -0.589119374752045,
    -0.59313553571701,
    0.398999333381653,
    -0.669425129890442,
    -0.0824997052550316,
    1.10110199451447,
    1.09411704540253,
    0.496543049812317,
    0.864346265792847,
    -0.000329156988300383,
    -1.39900064468384,
    -0.871215343475342,
    0.848287403583527,
    -0.450541794300079,
    -0.338643670082092,
    -0.230034217238426,
    -0.237949028611183,
    0.26269268989563,
    -0.315437704324722,
    -0.253818035125732,
    1.08865928649902,
    -0.076394334435463,
    0.566563367843628,
    0.00103327655233443,
    0.109523482620716,
    -0.468034714460373,
    0.663373827934265,
    1.03212201595306,
    -0.965081572532654,
    -0.355984181165695,
    0.481167733669281,
    0.939341366291046,
    -1.43014979362488,
    -0.308001816272736,
    0.500514805316925,
    -1.01503050327301,
    0.123037494719028,
    0.157860547304153,
    -1.20148706436157,
    -0.304913431406021,
    0.599206209182739,
    0.169961929321289,
    0.132926970720291,
    0.327926248311996,
    0.809473216533661,
    -0.842033326625824,
    -0.943092823028564,
    -0.0711575224995613,
    0.370576560497284,
    0.940639555454254,
    -0.25031653046608,
    0.601141631603241,
    -0.548416614532471,
    -0.305320858955383,
    0.28328338265419,
    0.239376679062843,
    0.193948447704315,
    -0.395745307207108,
    -0.987083733081818,
    0.818529903888702,
    -0.00782305188477039,
    -0.419385641813278,
    0.0678071975708008,
    -0.849334716796875,
    0.414339423179626,
    -0.874861598014832,
    0.350040763616562,
    0.125092908740044,
    -2.28567171096802,
    0.178147941827774,
    -1.24453294277191,
    -0.536106646060944,
    -0.33182767033577,
    0.216807708144188,
    0.962031126022339,
    0.266070246696472,
    0.353542536497116,
    -0.746428847312927,
    -0.112563341856003,
    0.431069135665894,
    0.210946083068848,
    0.136623814702034,
    -0.317520797252655,
    0.904275715351105,
    1.36134874820709,
    0.722247779369354,
    0.181408777832985,
    1.10303354263306,
    0.672251045703888,
    0.698297321796417,
    1.40209496021271,
    -0.364504516124725,
    0.0594567619264126,
    -0.135279625654221,
    -0.23151583969593,
    0.261116981506348,
    0.158680856227875,
    -0.550595581531525,
    0.129326865077019,
    0.498884439468384,
    0.775023221969604,
    -1.15530705451965,
    0.251861751079559,
    -0.824566602706909,
    -1.20348429679871,
    -0.361452966928482,
    1.0070903301239,
    -0.324127078056335,
    -0.254802525043488,
    0.240294933319092,
    -0.123521827161312,
    -1.06303858757019,
    0.624724924564362,
    -0.679011106491089,
    0.447645723819733,
    -0.398750871419907,
    -0.070178896188736,
    -0.921346306800842,
    0.473232746124268,
    -0.72577041387558,
    -1.43326234817505,
    0.325914084911346,
    0.572706580162048,
    -0.35582959651947,
    -0.940600395202637,
    -0.304235935211182,
    1.19818484783173,
    -2.14569902420044,
    -0.258589446544647,
    -1.05650389194489,
    0.245821982622147,
    -0.455048173666,
    0.111858673393726,
    -0.666933059692383,
    0.480834454298019,
    0.335937559604645,
    -0.593557775020599,
    -0.31837671995163,
    -0.597215116024017,
    0.334227859973907,
    0.321096271276474,
    -0.0689573511481285,
    0.433640897274017,
    -0.044014647603035,
    0.16673955321312,
    0.213584989309311,
    0.455057501792908,
    0.384808123111725,
    1.45602095127106,
    -1.76953423023224,
    -0.788201570510864,
    -1.39887094497681,
    -0.987034201622009,
    -0.308007955551147,
    0.0303755272179842,
    -1.13289976119995,
    -1.18198478221893,
    -0.16719862818718,
    0.867197036743164,
    -0.273155570030212,
    -1.44690632820129,
    -0.162217378616333,
    0.307093858718872,
    -0.526483416557312,
    -0.922611653804779,
    0.250266164541245,
    -0.048527754843235,
    -0.325054556131363,
    -0.295904666185379,
    0.401591271162033,
    2.37485480308533,
    0.671164512634277,
    -0.517287075519562,
    0.977664172649384,
    0.887623727321625,
    0.385577201843262,
    0.373436331748962,
    0.115243643522263,
    -0.405085206031799,
    0.617261648178101,
    -0.842770636081696,
    1.29267609119415,
    -0.368617177009583,
    0.987100780010223,
    1.18235564231873,
    0.127058207988739,
    -0.248720690608025,
    0.702932596206665,
    0.736054003238678,
    -0.151821777224541,
    0.16470967233181,
    1.43320500850677,
    0.803732216358185,
    -0.186717495322227,
    0.13795517385006,
    0.225798100233078 



// 0.928747475147247,
// -0.22091244161129,
// 0.340678840875626,
// 0.232478484511375,
// -0.0824997052550316,
// 1.08865928649902,
// 0.123037494719028,
// -0.548416614532471,
// -2.28567171096802,
// 1.36134874820709,
// 0.775023221969604,
// -0.070178896188736,
// 0.111858673393726,
// 1.45602095127106,
// -0.922611653804779,
// -0.842770636081696,
// 0.118580423295498,
// -0.443158328533173,
// 1.26924169063568,
// 0.27703258395195,
// 1.10110199451447,
// -0.076394334435463,
// 0.157860547304153,
// -0.305320858955383,
// 0.178147941827774,
// 0.722247779369354,
// -1.15530705451965,
// -0.921346306800842,
// -0.666933059692383,
// -1.76953423023224,
// 0.250266164541245,
// 1.29267609119415,
// 0.0733854696154594,
// -0.997133553028107,
// 0.763332605361938,
// 0.256476700305939,
// 1.09411704540253,
// 0.566563367843628,
// -1.20148706436157,
// 0.28328338265419,
// -1.24453294277191,
// 0.181408777832985,
// 0.251861751079559,
// 0.473232746124268,
// 0.480834454298019,
// -0.788201570510864,
// -0.048527754843235,
// -0.368617177009583,
// -0.216320112347603,
// -0.994887948036194,
// 1.23038351535797,
// -0.456453740596771,
// 0.496543049812317,
// 0.00103327655233443,
// -0.304913431406021,
// 0.239376679062843,
// -0.536106646060944,
// 1.10303354263306,
// -0.824566602706909,
// -0.72577041387558,
// 0.335937559604645,
// -1.39887094497681,
// -0.325054556131363,
// 0.987100780010223,
// 0.169354125857353,
// -0.131350949406624,
// 0.810515403747559,
// 0.390097618103027,
// 0.864346265792847,
// 0.109523482620716,
// 0.599206209182739,
// 0.193948447704315,
// -0.33182767033577,
// 0.672251045703888,
// -1.20348429679871,
// -1.43326234817505,
// -0.593557775020599,
// -0.987034201622009,
// -0.295904666185379,
// 1.18235564231873,
// -0.428392708301544,
// -0.341392010450363,
// 0.351482480764389,
// -0.0704145431518555,
// -0.000329156988300383,
// -0.468034714460373,
// 0.169961929321289,
// -0.395745307207108,
// 0.216807708144188,
// 0.698297321796417,
// -0.361452966928482,
// 0.325914084911346,
// -0.31837671995163,
// -0.308007955551147,
// 0.401591271162033,
// 0.127058207988739,
// 0.843611776828766,
// 0.456412494182587,
// 0.302027016878128,
// -1.4993953704834,
// -1.39900064468384,
// 0.663373827934265,
// 0.132926970720291,
// -0.987083733081818,
// 0.962031126022339,
// 1.40209496021271,
// 1.0070903301239,
// 0.572706580162048,
// -0.597215116024017,
// 0.0303755272179842,
// 2.37485480308533,
// -0.248720690608025,
// 0.00536061264574528,
// 0.00402470072731376,
// -0.450287461280823,
// -0.247053876519203,
// -0.871215343475342,
// 1.03212201595306,
// 0.327926248311996,
// 0.818529903888702,
// 0.266070246696472,
// -0.364504516124725,
// -0.324127078056335,
// -0.35582959651947,
// 0.334227859973907,
// -1.13289976119995,
// 0.671164512634277,
// 0.702932596206665,
// -0.137300744652748,
// 0.350713819265366,
// 1.04432547092438,
// -0.0885710567235947,
// 0.848287403583527,
// -0.965081572532654,
// 0.809473216533661,
// -0.00782305188477039,
// 0.353542536497116,
// 0.0594567619264126,
// -0.254802525043488,
// -0.940600395202637,
// 0.321096271276474,
// -1.18198478221893,
// -0.517287075519562,
// 0.736054003238678,
// -0.49064639210701,
// 0.262174665927887,
// 0.673821926116943,
// 1.01270878314972,
// -0.450541794300079,
// -0.355984181165695,
// -0.842033326625824,
// -0.419385641813278,
// -0.746428847312927,
// -0.135279625654221,
// 0.240294933319092,
// -0.304235935211182,
// -0.0689573511481285,
// -0.16719862818718,
// 0.977664172649384,
// -0.151821777224541,
// 0.230988517403603,
// -0.048433355987072,
// -0.16774320602417,
// -0.245661944150925,
// -0.338643670082092,
// 0.481167733669281,
// -0.943092823028564,
// 0.0678071975708008,
// -0.112563341856003,
// -0.23151583969593,
// -0.123521827161312,
// 1.19818484783173,
// 0.433640897274017,
// 0.867197036743164,
// 0.887623727321625,
// 0.16470967233181,
// -1.13217556476593,
// 0.543828129768372,
// 1.87203669548035,
// -0.478165984153748,
// -0.230034217238426,
// 0.939341366291046,
// -0.0711575224995613,
// -0.849334716796875,
// 0.431069135665894,
// 0.261116981506348,
// -1.06303858757019,
// -2.14569902420044,
// -0.044014647603035,
// -0.273155570030212,
// 0.385577201843262,
// 1.43320500850677,
// 0.981262445449829,
// 0.915902376174927,
// 0.0531666539609432,
// -0.589119374752045,
// -0.237949028611183,
// -1.43014979362488,
// 0.370576560497284,
// 0.414339423179626,
// 0.210946083068848,
// 0.158680856227875,
// 0.624724924564362,
// -0.258589446544647,
// 0.16673955321312,
// -1.44690632820129,
// 0.373436331748962,
// 0.803732216358185,
// -0.275694042444229,
// 0.681087553501129,
// 0.374954998493195,
// -0.59313553571701,
// 0.26269268989563,
// -0.308001816272736,
// 0.940639555454254,
// -0.874861598014832,
// 0.136623814702034,
// -0.550595581531525,
// -0.679011106491089,
// -1.05650389194489,
// 0.213584989309311,
// -0.162217378616333,
// 0.115243643522263,
// -0.186717495322227,
// 0.458381682634354,
// -0.425360172986984,
// 0.566317498683929,
// 0.398999333381653,
// -0.315437704324722,
// 0.500514805316925,
// -0.25031653046608,
// 0.350040763616562,
// -0.317520797252655,
// 0.129326865077019,
// 0.447645723819733,
// 0.245821982622147,
// 0.455057501792908,
// 0.307093858718872,
// -0.405085206031799,
// 0.13795517385006,
// -0.859836757183075,
// 0.208959579467773,
// 0.792478382587433,
// -0.669425129890442,
// -0.253818035125732,
// -1.01503050327301,
// 0.601141631603241,
// 0.125092908740044,
// 0.904275715351105,
// 0.498884439468384,
// -0.398750871419907,
// -0.455048173666,
// 0.384808123111725,
// -0.526483416557312,
// 0.617261648178101,
// 0.225798100233078
};

constexpr std::array<double, 264> W4_data =
{
    1.15496790409088,
    0.720434963703156,
    1.07307076454163,
    -2.68247699737549,
    -2.64888548851013,
    0.259213000535965,
    1.51967513561249,
    1.54009389877319,
    -1.04019522666931,
    2.10143280029297,
    0.0043989447876811,
    -1.25469100475311,
    -2.51979541778564,
    0.901214301586151,
    0.847959101200104,
    1.30571949481964,
    -2.17798280715942,
    -0.9259153008461,
    1.84006333351135,
    0.789070963859558,
    -0.370641678571701,
    0.194479271769524,
    -0.378281027078629,
    0.925322592258453,
    1.26704525947571,
    -1.09689104557037,
    2.09795904159546,
    -2.0838987827301,
    0.11978618055582,
    -1.21036434173584,
    1.68464481830597,
    -0.135086610913277,
    -0.999509394168854,
    -1.60115313529968,
    0.934219360351563,
    1.28817009925842,
    -0.916245460510254,
    1.27218496799469,
    -0.998437225818634,
    0.471503436565399,
    1.34095215797424,
    0.860655426979065,
    -0.760026633739471,
    -0.800763785839081,
    -0.649541139602661,
    -0.343869388103485,
    -1.70810973644257,
    2.4398946762085,
    0.0827361866831779,
    -0.490477532148361,
    -1.34021723270416,
    2.03310942649841,
    0.0823729187250137,
    -1.46317648887634,
    0.738475382328033,
    1.25785028934479,
    0.459993958473206,
    0.929783463478088,
    -0.0812694504857063,
    -0.640182197093964,
    0.202272862195969,
    -1.12464797496796,
    -0.307410687208176,
    1.79511415958405



// 1.15496790409088,
// -2.64888548851013,
// -1.04019522666931,
// -2.51979541778564,
// -2.17798280715942,
// -0.370641678571701,
// 1.26704525947571,
// 0.11978618055582,
// -0.999509394168854,
// -0.916245460510254,
// 1.34095215797424,
// -0.649541139602661,
// 0.0827361866831779,
// 0.0823729187250137,
// 0.459993958473206,
// 0.202272862195969,
// 0.720434963703156,
// 0.259213000535965,
// 2.10143280029297,
// 0.901214301586151,
// -0.9259153008461,
// 0.194479271769524,
// -1.09689104557037,
// -1.21036434173584,
// -1.60115313529968,
// 1.27218496799469,
// 0.860655426979065,
// -0.343869388103485,
// -0.490477532148361,
// -1.46317648887634,
// 0.929783463478088,
// -1.12464797496796,
// 1.07307076454163,
// 1.51967513561249,
// 0.0043989447876811,
// 0.847959101200104,
// 1.84006333351135,
// -0.378281027078629,
// 2.09795904159546,
// 1.68464481830597,
// 0.934219360351563,
// -0.998437225818634,
// -0.760026633739471,
// -1.70810973644257,
// -1.34021723270416,
// 0.738475382328033,
// -0.0812694504857063,
// -0.307410687208176,
// -2.68247699737549,
// 1.54009389877319,
// -1.25469100475311,
// 1.30571949481964,
// 0.789070963859558,
// 0.925322592258453,
// -2.0838987827301,
// -0.135086610913277,
// 1.28817009925842,
// 0.471503436565399,
// -0.800763785839081,
// 2.4398946762085,
// 2.03310942649841,
// 1.25785028934479,
// -0.640182197093964,
// 1.79511415958405
};

constexpr std::array<double, 16> B1_data =
{
    -0.856053233146667,
    0.299127906560898,
    -0.520716726779938,
    0.772461295127869,
    0.0834392830729485,
    -0.405647426843643,
    0.343869537115097,
    -0.848781526088715,
    1.04696619510651,
    0.775198459625244,
    -0.184610620141029,
    -0.597316205501556,
    -0.132456228137016,
    0.0218230783939362,
    1.55549299716949,
    0.582780420780182
};

constexpr std::array<double, 16> B2_data =
{
    0.40303909778595,
    -0.793079137802124,
    0.221006974577904,
    0.158226385712624,
    0.39412197470665,
    0.382671117782593,
    -0.1337820738554,
    0.0864124000072479,
    -0.455804467201233,
    0.712578654289246,
    0.372920542955399,
    0.847211480140686,
    -0.267316162586212,
    0.632389307022095,
    -0.109682455658913,
    -0.142096906900406
};

constexpr std::array<double, 16> B3_data =
{
    0.853201389312744,
    0.32586207985878,
    0.540267407894135,
    0.288121372461319,
    0.0801954492926598,
    -0.147755175828934,
    0.352727681398392,
    -0.309247881174088,
    -0.0521119013428688,
    -0.00933279749006033,
    0.287271529436111,
    -0.326762765645981,
    -0.442708432674408,
    -0.50548255443573,
    0.310500890016556,
    -0.0361274145543575
};

constexpr std::array<double, 4> B4_data =
{
    0.265564352273941,
    0.342342644929886,
    0.0911224335432053,
    -0.699481189250946    
};


#endif