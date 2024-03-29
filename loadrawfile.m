function tfringeMat = loadrawfile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   TFRINGEMAT = IMPORTFILE(FILENAME)
%   Reads data from text file FILENAME for the default selection.
%
%   TFRINGEMAT = IMPORTFILE(FILENAME, STARTROW, ENDROW)
%   Reads data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   tfringeMat = importfile('CBI_AA Merge.raw.txt', 3, 9302);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2022/05/12 16:07:24

%% Initialize variables.
if nargin<=2
    startRow = 3;
    endRow = inf;
end

%% Format for each line of text:
%   column6: double (%f)
%	column7: double (%f)
%   column8: double (%f)
%	column9: double (%f)
%   column10: double (%f)
%	column11: double (%f)
%   column12: double (%f)
%	column13: double (%f)
%   column14: double (%f)
%	column15: double (%f)
%   column16: double (%f)
%	column17: double (%f)
%   column18: double (%f)
%	column19: double (%f)
%   column20: double (%f)
%	column21: double (%f)
%   column22: double (%f)
%	column23: double (%f)
%   column24: double (%f)
%	column25: double (%f)
%   column26: double (%f)
%	column27: double (%f)
%   column28: double (%f)
%	column29: double (%f)
%   column30: double (%f)
%	column31: double (%f)
%   column32: double (%f)
%	column33: double (%f)
%   column34: double (%f)
%	column35: double (%f)
%   column36: double (%f)
%	column37: double (%f)
%   column38: double (%f)
%	column39: double (%f)
%   column40: double (%f)
%	column41: double (%f)
%   column42: double (%f)
%	column43: double (%f)
%   column44: double (%f)
%	column45: double (%f)
%   column46: double (%f)
%	column47: double (%f)
%   column48: double (%f)
%	column49: double (%f)
%   column50: double (%f)
%	column51: double (%f)
%   column52: double (%f)
%	column53: double (%f)
%   column54: double (%f)
%	column55: double (%f)
%   column56: double (%f)
%	column57: double (%f)
%   column58: double (%f)
%	column59: double (%f)
%   column60: double (%f)
%	column61: double (%f)
%   column62: double (%f)
%	column63: double (%f)
%   column64: double (%f)
%	column65: double (%f)
%   column66: double (%f)
%	column67: double (%f)
%   column68: double (%f)
%	column69: double (%f)
%   column70: double (%f)
%	column71: double (%f)
%   column72: double (%f)
%	column73: double (%f)
%   column74: double (%f)
%	column75: double (%f)
%   column76: double (%f)
%	column77: double (%f)
%   column78: double (%f)
%	column79: double (%f)
%   column80: double (%f)
%	column81: double (%f)
%   column82: double (%f)
%	column83: double (%f)
%   column84: double (%f)
%	column85: double (%f)
%   column86: double (%f)
%	column87: double (%f)
%   column88: double (%f)
%	column89: double (%f)
%   column90: double (%f)
%	column91: double (%f)
%   column92: double (%f)
%	column93: double (%f)
%   column94: double (%f)
%	column95: double (%f)
%   column96: double (%f)
%	column97: double (%f)
%   column98: double (%f)
%	column99: double (%f)
%   column100: double (%f)
%	column101: double (%f)
%   column102: double (%f)
%	column103: double (%f)
%   column104: double (%f)
%	column105: double (%f)
%   column106: double (%f)
%	column107: double (%f)
%   column108: double (%f)
%	column109: double (%f)
%   column110: double (%f)
%	column111: double (%f)
%   column112: double (%f)
%	column113: double (%f)
%   column114: double (%f)
%	column115: double (%f)
%   column116: double (%f)
%	column117: double (%f)
%   column118: double (%f)
%	column119: double (%f)
%   column120: double (%f)
%	column121: double (%f)
%   column122: double (%f)
%	column123: double (%f)
%   column124: double (%f)
%	column125: double (%f)
%   column126: double (%f)
%	column127: double (%f)
%   column128: double (%f)
%	column129: double (%f)
%   column130: double (%f)
%	column131: double (%f)
%   column132: double (%f)
%	column133: double (%f)
%   column134: double (%f)
%	column135: double (%f)
%   column136: double (%f)
%	column137: double (%f)
%   column138: double (%f)
%	column139: double (%f)
%   column140: double (%f)
%	column141: double (%f)
%   column142: double (%f)
%	column143: double (%f)
%   column144: double (%f)
%	column145: double (%f)
%   column146: double (%f)
%	column147: double (%f)
%   column148: double (%f)
%	column149: double (%f)
%   column150: double (%f)
%	column151: double (%f)
%   column152: double (%f)
%	column153: double (%f)
%   column154: double (%f)
%	column155: double (%f)
%   column156: double (%f)
%	column157: double (%f)
%   column158: double (%f)
%	column159: double (%f)
%   column160: double (%f)
%	column161: double (%f)
%   column162: double (%f)
%	column163: double (%f)
%   column164: double (%f)
%	column165: double (%f)
%   column166: double (%f)
%	column167: double (%f)
%   column168: double (%f)
%	column169: double (%f)
%   column170: double (%f)
%	column171: double (%f)
%   column172: double (%f)
%	column173: double (%f)
%   column174: double (%f)
%	column175: double (%f)
%   column176: double (%f)
%	column177: double (%f)
%   column178: double (%f)
%	column179: double (%f)
%   column180: double (%f)
%	column181: double (%f)
%   column182: double (%f)
%	column183: double (%f)
%   column184: double (%f)
%	column185: double (%f)
%   column186: double (%f)
%	column187: double (%f)
%   column188: double (%f)
%	column189: double (%f)
%   column190: double (%f)
%	column191: double (%f)
%   column192: double (%f)
%	column193: double (%f)
%   column194: double (%f)
%	column195: double (%f)
%   column196: double (%f)
%	column197: double (%f)
%   column198: double (%f)
%	column199: double (%f)
%   column200: double (%f)
%	column201: double (%f)
%   column202: double (%f)
%	column203: double (%f)
%   column204: double (%f)
%	column205: double (%f)
%   column206: double (%f)
%	column207: double (%f)
%   column208: double (%f)
%	column209: double (%f)
%   column210: double (%f)
%	column211: double (%f)
%   column212: double (%f)
%	column213: double (%f)
%   column214: double (%f)
%	column215: double (%f)
%   column216: double (%f)
%	column217: double (%f)
%   column218: double (%f)
%	column219: double (%f)
%   column220: double (%f)
%	column221: double (%f)
%   column222: double (%f)
%	column223: double (%f)
%   column224: double (%f)
%	column225: double (%f)
%   column226: double (%f)
%	column227: double (%f)
%   column228: double (%f)
%	column229: double (%f)
%   column230: double (%f)
%	column231: double (%f)
%   column232: double (%f)
%	column233: double (%f)
%   column234: double (%f)
%	column235: double (%f)
%   column236: double (%f)
%	column237: double (%f)
%   column238: double (%f)
%	column239: double (%f)
%   column240: double (%f)
%	column241: double (%f)
%   column242: double (%f)
%	column243: double (%f)
%   column244: double (%f)
%	column245: double (%f)
%   column246: double (%f)
%	column247: double (%f)
%   column248: double (%f)
%	column249: double (%f)
%   column250: double (%f)
%	column251: double (%f)
%   column252: double (%f)
%	column253: double (%f)
%   column254: double (%f)
%	column255: double (%f)
%   column256: double (%f)
%	column257: double (%f)
%   column258: double (%f)
%	column259: double (%f)
%   column260: double (%f)
%	column261: double (%f)
%   column262: double (%f)
%	column263: double (%f)
%   column264: double (%f)
%	column265: double (%f)
%   column266: double (%f)
%	column267: double (%f)
%   column268: double (%f)
%	column269: double (%f)
%   column270: double (%f)
%	column271: double (%f)
%   column272: double (%f)
%	column273: double (%f)
%   column274: double (%f)
%	column275: double (%f)
%   column276: double (%f)
%	column277: double (%f)
%   column278: double (%f)
%	column279: double (%f)
%   column280: double (%f)
%	column281: double (%f)
%   column282: double (%f)
%	column283: double (%f)
%   column284: double (%f)
%	column285: double (%f)
%   column286: double (%f)
%	column287: double (%f)
%   column288: double (%f)
%	column289: double (%f)
%   column290: double (%f)
%	column291: double (%f)
%   column292: double (%f)
%	column293: double (%f)
%   column294: double (%f)
%	column295: double (%f)
%   column296: double (%f)
%	column297: double (%f)
%   column298: double (%f)
%	column299: double (%f)
%   column300: double (%f)
%	column301: double (%f)
%   column302: double (%f)
%	column303: double (%f)
%   column304: double (%f)
%	column305: double (%f)
%   column306: double (%f)
%	column307: double (%f)
%   column308: double (%f)
%	column309: double (%f)
%   column310: double (%f)
%	column311: double (%f)
%   column312: double (%f)
%	column313: double (%f)
%   column314: double (%f)
%	column315: double (%f)
%   column316: double (%f)
%	column317: double (%f)
%   column318: double (%f)
%	column319: double (%f)
%   column320: double (%f)
%	column321: double (%f)
%   column322: double (%f)
%	column323: double (%f)
%   column324: double (%f)
%	column325: double (%f)
%   column326: double (%f)
%	column327: double (%f)
%   column328: double (%f)
%	column329: double (%f)
%   column330: double (%f)
%	column331: double (%f)
%   column332: double (%f)
%	column333: double (%f)
%   column334: double (%f)
%	column335: double (%f)
%   column336: double (%f)
%	column337: double (%f)
%   column338: double (%f)
%	column339: double (%f)
%   column340: double (%f)
%	column341: double (%f)
%   column342: double (%f)
%	column343: double (%f)
%   column344: double (%f)
%	column345: double (%f)
%   column346: double (%f)
%	column347: double (%f)
%   column348: double (%f)
%	column349: double (%f)
%   column350: double (%f)
%	column351: double (%f)
%   column352: double (%f)
%	column353: double (%f)
%   column354: double (%f)
%	column355: double (%f)
%   column356: double (%f)
%	column357: double (%f)
%   column358: double (%f)
%	column359: double (%f)
%   column360: double (%f)
%	column361: double (%f)
%   column362: double (%f)
%	column363: double (%f)
%   column364: double (%f)
%	column365: double (%f)
%   column366: double (%f)
%	column367: double (%f)
%   column368: double (%f)
%	column369: double (%f)
%   column370: double (%f)
%	column371: double (%f)
%   column372: double (%f)
%	column373: double (%f)
%   column374: double (%f)
%	column375: double (%f)
%   column376: double (%f)
%	column377: double (%f)
%   column378: double (%f)
%	column379: double (%f)
%   column380: double (%f)
%	column381: double (%f)
%   column382: double (%f)
%	column383: double (%f)
%   column384: double (%f)
%	column385: double (%f)
%   column386: double (%f)
%	column387: double (%f)
%   column388: double (%f)
%	column389: double (%f)
%   column390: double (%f)
%	column391: double (%f)
%   column392: double (%f)
%	column393: double (%f)
%   column394: double (%f)
%	column395: double (%f)
%   column396: double (%f)
%	column397: double (%f)
%   column398: double (%f)
%	column399: double (%f)
%   column400: double (%f)
%	column401: double (%f)
%   column402: double (%f)
%	column403: double (%f)
%   column404: double (%f)
%	column405: double (%f)
%   column406: double (%f)
%	column407: double (%f)
%   column408: double (%f)
%	column409: double (%f)
%   column410: double (%f)
%	column411: double (%f)
%   column412: double (%f)
%	column413: double (%f)
%   column414: double (%f)
%	column415: double (%f)
%   column416: double (%f)
%	column417: double (%f)
%   column418: double (%f)
%	column419: double (%f)
%   column420: double (%f)
%	column421: double (%f)
%   column422: double (%f)
%	column423: double (%f)
%   column424: double (%f)
%	column425: double (%f)
%   column426: double (%f)
%	column427: double (%f)
%   column428: double (%f)
%	column429: double (%f)
%   column430: double (%f)
%	column431: double (%f)
%   column432: double (%f)
%	column433: double (%f)
%   column434: double (%f)
%	column435: double (%f)
%   column436: double (%f)
%	column437: double (%f)
%   column438: double (%f)
%	column439: double (%f)
%   column440: double (%f)
%	column441: double (%f)
%   column442: double (%f)
%	column443: double (%f)
%   column444: double (%f)
%	column445: double (%f)
%   column446: double (%f)
%	column447: double (%f)
%   column448: double (%f)
%	column449: double (%f)
%   column450: double (%f)
%	column451: double (%f)
%   column452: double (%f)
%	column453: double (%f)
%   column454: double (%f)
%	column455: double (%f)
%   column456: double (%f)
%	column457: double (%f)
%   column458: double (%f)
%	column459: double (%f)
%   column460: double (%f)
%	column461: double (%f)
%   column462: double (%f)
%	column463: double (%f)
%   column464: double (%f)
%	column465: double (%f)
%   column466: double (%f)
%	column467: double (%f)
%   column468: double (%f)
%	column469: double (%f)
%   column470: double (%f)
%	column471: double (%f)
%   column472: double (%f)
%	column473: double (%f)
%   column474: double (%f)
%	column475: double (%f)
%   column476: double (%f)
%	column477: double (%f)
%   column478: double (%f)
%	column479: double (%f)
%   column480: double (%f)
%	column481: double (%f)
%   column482: double (%f)
%	column483: double (%f)
%   column484: double (%f)
%	column485: double (%f)
%   column486: double (%f)
%	column487: double (%f)
%   column488: double (%f)
%	column489: double (%f)
%   column490: double (%f)
%	column491: double (%f)
%   column492: double (%f)
%	column493: double (%f)
%   column494: double (%f)
%	column495: double (%f)
%   column496: double (%f)
%	column497: double (%f)
%   column498: double (%f)
%	column499: double (%f)
%   column500: double (%f)
%	column501: double (%f)
%   column502: double (%f)
%	column503: double (%f)
%   column504: double (%f)
%	column505: double (%f)
%   column506: double (%f)
%	column507: double (%f)
%   column508: double (%f)
%	column509: double (%f)
%   column510: double (%f)
%	column511: double (%f)
%   column512: double (%f)
%	column513: double (%f)
%   column514: double (%f)
%	column515: double (%f)
%   column516: double (%f)
%	column517: double (%f)
%   column518: double (%f)
%	column519: double (%f)
%   column520: double (%f)
%	column521: double (%f)
%   column522: double (%f)
%	column523: double (%f)
%   column524: double (%f)
%	column525: double (%f)
%   column526: double (%f)
%	column527: double (%f)
%   column528: double (%f)
%	column529: double (%f)
%   column530: double (%f)
%	column531: double (%f)
%   column532: double (%f)
%	column533: double (%f)
%   column534: double (%f)
%	column535: double (%f)
%   column536: double (%f)
%	column537: double (%f)
%   column538: double (%f)
%	column539: double (%f)
%   column540: double (%f)
%	column541: double (%f)
%   column542: double (%f)
%	column543: double (%f)
%   column544: double (%f)
%	column545: double (%f)
%   column546: double (%f)
%	column547: double (%f)
%   column548: double (%f)
%	column549: double (%f)
%   column550: double (%f)
%	column551: double (%f)
%   column552: double (%f)
%	column553: double (%f)
%   column554: double (%f)
%	column555: double (%f)
%   column556: double (%f)
%	column557: double (%f)
%   column558: double (%f)
%	column559: double (%f)
%   column560: double (%f)
%	column561: double (%f)
%   column562: double (%f)
%	column563: double (%f)
%   column564: double (%f)
%	column565: double (%f)
%   column566: double (%f)
%	column567: double (%f)
%   column568: double (%f)
%	column569: double (%f)
%   column570: double (%f)
%	column571: double (%f)
%   column572: double (%f)
%	column573: double (%f)
%   column574: double (%f)
%	column575: double (%f)
%   column576: double (%f)
%	column577: double (%f)
%   column578: double (%f)
%	column579: double (%f)
%   column580: double (%f)
%	column581: double (%f)
%   column582: double (%f)
%	column583: double (%f)
%   column584: double (%f)
%	column585: double (%f)
%   column586: double (%f)
%	column587: double (%f)
%   column588: double (%f)
%	column589: double (%f)
%   column590: double (%f)
%	column591: double (%f)
%   column592: double (%f)
%	column593: double (%f)
%   column594: double (%f)
%	column595: double (%f)
%   column596: double (%f)
%	column597: double (%f)
%   column598: double (%f)
%	column599: double (%f)
%   column600: double (%f)
%	column601: double (%f)
%   column602: double (%f)
%	column603: double (%f)
%   column604: double (%f)
%	column605: double (%f)
%   column606: double (%f)
%	column607: double (%f)
%   column608: double (%f)
%	column609: double (%f)
%   column610: double (%f)
%	column611: double (%f)
%   column612: double (%f)
%	column613: double (%f)
%   column614: double (%f)
%	column615: double (%f)
%   column616: double (%f)
%	column617: double (%f)
%   column618: double (%f)
%	column619: double (%f)
%   column620: double (%f)
%	column621: double (%f)
%   column622: double (%f)
%	column623: double (%f)
%   column624: double (%f)
%	column625: double (%f)
%   column626: double (%f)
%	column627: double (%f)
%   column628: double (%f)
%	column629: double (%f)
%   column630: double (%f)
%	column631: double (%f)
%   column632: double (%f)
%	column633: double (%f)
%   column634: double (%f)
%	column635: double (%f)
%   column636: double (%f)
%	column637: double (%f)
%   column638: double (%f)
%	column639: double (%f)
%   column640: double (%f)
%	column641: double (%f)
%   column642: double (%f)
%	column643: double (%f)
%   column644: double (%f)
%	column645: double (%f)
%   column646: double (%f)
%	column647: double (%f)
%   column648: double (%f)
%	column649: double (%f)
%   column650: double (%f)
%	column651: double (%f)
%   column652: double (%f)
%	column653: double (%f)
%   column654: double (%f)
%	column655: double (%f)
%   column656: double (%f)
%	column657: double (%f)
%   column658: double (%f)
%	column659: double (%f)
%   column660: double (%f)
%	column661: double (%f)
%   column662: double (%f)
%	column663: double (%f)
%   column664: double (%f)
%	column665: double (%f)
%   column666: double (%f)
%	column667: double (%f)
%   column668: double (%f)
%	column669: double (%f)
%   column670: double (%f)
%	column671: double (%f)
%   column672: double (%f)
%	column673: double (%f)
%   column674: double (%f)
%	column675: double (%f)
%   column676: double (%f)
%	column677: double (%f)
%   column678: double (%f)
%	column679: double (%f)
%   column680: double (%f)
%	column681: double (%f)
%   column682: double (%f)
%	column683: double (%f)
%   column684: double (%f)
%	column685: double (%f)
%   column686: double (%f)
%	column687: double (%f)
%   column688: double (%f)
%	column689: double (%f)
%   column690: double (%f)
%	column691: double (%f)
%   column692: double (%f)
%	column693: double (%f)
%   column694: double (%f)
%	column695: double (%f)
%   column696: double (%f)
%	column697: double (%f)
%   column698: double (%f)
%	column699: double (%f)
%   column700: double (%f)
%	column701: double (%f)
%   column702: double (%f)
%	column703: double (%f)
%   column704: double (%f)
%	column705: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*25s%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%15f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post processing code is included. To generate code which works for unimportable data, select unimportable cells in a file and regenerate the script.

%% Create output variable
tfringeMat = [dataArray{1:end-1}];
