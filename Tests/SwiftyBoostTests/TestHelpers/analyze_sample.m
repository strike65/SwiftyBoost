
function out = analyze_sample(csvA, varargin)
% ANALYZE_SAMPLE  Kennzahlen aus Stichprobe(n) in CSV berechnen.
%
% Aufruf (stetig, 1D):
%   out = analyze_sample('sample.csv', 'type','continuous');
%
% Aufruf (diskret):
%   out = analyze_sample('counts.csv', 'type','discrete');
%
% Aufruf mit zweiter Stichprobe für KL-Divergenz:
%   out = analyze_sample('p.csv', 'type','continuous', 'csvB','q.csv');
%
% Wichtige Optionen (Name-Wert-Paare):
%   'type'        : 'discrete' | 'continuous'  (Default: automatische Heuristik)
%   'column'      : Spaltenindex für 1D-Auswertung (Default: 1)
%   'entropyBase' : 2 (Bits) | exp(1) (Nats) | 10 (Hartleys) (Default: exp(1))
%   'bins'        : Anzahl Bins (diskret) oder für Histogrammodusschätzer (Default: auto)
%   'bandwidth'   : KDE-Bandbreite h (wenn leer -> Silverman) (nur stetig)
%   'k'           : k für kNN-Entropie/KL (Default: 5)
%   'csvB'        : zweites CSV für KL-Divergenz P||Q (optional)
%   'millerMadow' : true/false Korrektur (diskrete Entropie) (Default: true)
%
% Rückgabe (Struktur):
%   out.mean, out.var, out.skewness, out.kurtosis_excess, out.median, out.mode
%   out.entropy, out.kl (falls csvB gesetzt), out.notes (Hinweise)
%
% Hinweise:
%   - Diskret: Entropie = Shannon (Plug-in), optional Miller-Madow.
%   - Stetig:  Differentialentropie via kNN (Kozachenko-Leonenko).
%   - Modus:   Diskret = häufigster Wert; stetig = Argmax KDE.
%   - KL:      Diskret per Häufigkeiten (gleicher Träger, mit Glättung);
%              stetig per kNN-Schätzer (1D/multivariat).
%
% Autor: (c) 2025

%% Optionen
p = inputParser;
addRequired(p, 'csvA', @(s)ischar(s) || isstring(s));
addParameter(p, 'type', '', @(s)ischar(s) || isstring(s));
addParameter(p, 'column', 1, @(x)isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'entropyBase', exp(1), @(x)isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'bins', [], @(x)isempty(x) || (isscalar(x) && x>=2));
addParameter(p, 'bandwidth', [], @(x)isempty(x) || (isscalar(x) && x>0));
addParameter(p, 'k', 5, @(x)isnumeric(x) && isscalar(x) && x>=1);
addParameter(p, 'csvB', '', @(s)ischar(s) || isstring(s));
addParameter(p, 'millerMadow', true, @(x)islogical(x));
parse(p, csvA, varargin{:});
opt = p.Results;

baseLog = @(x) log(x) / log(opt.entropyBase);

%% Daten laden
A = readmatrix(csvA);
A = A(~any(isnan(A),2),:);  % Zeilen mit NaN verwerfen
notes = {};

if isempty(A)
    error('Leere oder nur NaN-Zeilen in %s.', csvA);
end

% Typheuristik (falls nicht gesetzt)
if isempty(opt.type)
    % Wenn viele Duplikate und ganzzahlige Werte -> diskret
    col = min(opt.column, size(A,2));
    x0  = A(:,col);
    is_int = all(abs(x0 - round(x0)) < 1e-12);
    uniq_ratio = numel(unique(x0)) / numel(x0);
    if is_int && uniq_ratio < 0.25
        opt.type = 'discrete';
        notes{end+1} = 'Heuristik: Daten als diskret interpretiert.';
    else
        opt.type = 'continuous';
        notes{end+1} = 'Heuristik: Daten als stetig interpretiert.';
    end
else
    opt.type = lower(string(opt.type));
end

%% Ein-/mehrdimensional behandeln
d = size(A,2);

% Für skalare Lage-/Streuungsmaße nehmen wir standardmäßig eine Spalte
col = min(opt.column, d);
x  = A(:,col);

%% Grundmaße (pro Spalte, 1D)
out.mean  = mean(x);
out.var   = var(x, 1);           % Populationsvarianz; für Stichprobenvarianz var(x,0)
out.median = median(x);

% Bias-korrigierte Schiefe/Kurtosis (Exzess)
n = numel(x);
xm = x - mean(x);
s2 = sum(xm.^2) / (n-1);
s  = sqrt(s2);
m3 = sum(xm.^3) / n;
m4 = sum(xm.^4) / n;

if n>=3 && s>0
    g1 = (n / ((n-1)*(n-2))) * (m3 / (s^3));
else
    g1 = NaN;
    notes{end+1} = 'Schiefe nicht definiert (n<3 oder s=0).';
end

if n>=4 && s>0
    g2 = (n*(n+1)/((n-1)*(n-2)*(n-3))) * (m4 / (s^4)) - (3*(n-1)^2)/((n-2)*(n-3));
else
    g2 = NaN;
    notes{end+1} = 'Kurtosis nicht definiert (n<4 oder s=0).';
end

out.skewness = g1;
out.kurtosis_excess = g2;

%% Modus & Entropie je nach Typ
switch opt.type
    case 'discrete'
        % Diskrete Häufigkeiten auf Basis der ausgewählten Spalte
        [vals,~,idx] = unique(x);
        counts = accumarray(idx, 1);
        p = counts / sum(counts);
        % Modus: häufigster Wert
        [~,I] = max(p);
        out.mode = vals(I);

        % Entropie (Shannon, Plug-in), optional Miller–Madow
        H = -sum(p .* baseLog(p + (p==0))); % 0*log(0) = 0
        if opt.millerMadow
            K = sum(p>0);
            H = H + (K-1)/(2*n) / log(opt.entropyBase);
            notes{end+1} = 'Entropie: Plug-in + Miller–Madow-Korrektur.';
        else
            notes{end+1} = 'Entropie: Plug-in-Schätzer (ohne Bias-Korrektur).';
        end
        out.entropy = H;

    case 'continuous'
        % Modus via KDE (1D) auf der gewählten Spalte
        if isempty(opt.bandwidth)
            h = silverman_bandwidth(x);
            notes{end+1} = sprintf('KDE-Bandbreite (Silverman): h=%.6g', h);
        else
            h = opt.bandwidth;
        end
        gridN = max(512, min(4096, 2^nextpow2(n)));
        xs = linspace(min(x), max(x), gridN).';
        f  = kde_pdf_1d(x, xs, h);
        [~,I] = max(f);
        out.mode = xs(I);

        % Differentialentropie via Kozachenko–Leonenko (1D oder d>1)
        Hc = entropy_knn(A, opt.k); % nutzt alle Spalten (multivariat)
        out.entropy = Hc / log(opt.entropyBase);
        notes{end+1} = sprintf('Differentialentropie: kNN (k=%d). Basis=%g.', opt.k, opt.entropyBase);

    otherwise
        error('Unbekannter Typ: %s', opt.type);
end

%% KL-Divergenz (optional, P||Q)
out.kl = [];
if ~isempty(opt.csvB)
    B = readmatrix(opt.csvB);
    B = B(~any(isnan(B),2),:);
    if isempty(B)
        warning('Zweite CSV %s leer oder nur NaN-Zeilen – KL wird übersprungen.', opt.csvB);
    else
        switch opt.type
            case 'discrete'
                % Diskret: gleiche Trägermenge bilden, Laplace-Glättung
                y = B(:,col);
                [vals,~,idxP] = unique(x);
                countsP = accumarray(idxP, 1, [numel(vals),1]);
                % Zwinge Q auf denselben Träger
                [~, loc] = ismember(B(:,col), vals);
                countsQ = accumarray(max(loc,1), 1, [numel(vals),1]); % alle Nicht-Treffer auf 1. Index addieren -> wird geglättet
                alpha = 0.5;  % Add-alpha-Glättung
                pP = (countsP + alpha) / (sum(countsP) + alpha*numel(vals));
                pQ = (countsQ + alpha) / (sum(countsQ) + alpha*numel(vals));
                D = sum( pP .* baseLog(pP ./ pQ) );
                out.kl = D;
                notes{end+1} = sprintf('KL (diskret) mit Add-α (α=%.2g) und gemeinsamem Träger.', alpha);

            case 'continuous'
                % Stetig: kNN-KL-Schätzer (multivariat, gleicher d)
                if size(B,2) ~= size(A,2)
                    error('Dimensionalität P (%d) und Q (%d) differiert.', size(A,2), size(B,2));
                end
                D = kl_knn(A, B, opt.k);  % D_KL(P||Q)
                out.kl = D / log(opt.entropyBase);
                notes{end+1} = sprintf('KL (stetig) via kNN (k=%d). Basis=%g.', opt.k, opt.entropyBase);
        end
    end
end

out.notes = notes;

%% Ausgabe kurz auf Konsole
fprintf('n=%d, Spalte=%d, Typ=%s\n', n, col, string(opt.type));
fprintf('mean=%.6g, var=%.6g, median=%.6g, mode=%.6g\n', out.mean, out.var, out.median, out.mode);
fprintf('skew=%.6g, excess kurtosis=%.6g\n', out.skewness, out.kurtosis_excess);
fprintf('entropy=%.6g', out.entropy);
if ~isempty(out.kl), fprintf(', KL=%.6g', out.kl); end
fprintf('\n');

end % main function

%% ===== Hilfsfunktionen =====

function h = silverman_bandwidth(x)
% Silverman’s rule of thumb für 1D
x = x(:);
n = numel(x);
sigma = std(x, 0);
iqr_ = iqr(x);
s = min(sigma, iqr_/1.349);
h = 0.9 * s * n^(-1/5);
if h<=0 || ~isfinite(h)
    h = 1.06 * sigma * n^(-1/5); % Scott/Silverman fallback
end
end

function f = kde_pdf_1d(x, xs, h)
% Gaussian KDE (1D)
x  = x(:)';             % 1×n
xs = xs(:);             % m×1
u  = (xs - x) ./ h;     % m×n (implizit via bsxfun)
K  = exp(-0.5 * (u.^2)) / sqrt(2*pi);
f  = mean(K, 2) / h;    % m×1
end

function H = entropy_knn(X, k)
% Kozachenko–Leonenko-Differentialentropie für beliebige Dimension d
% Robust, ohne externe Helper; mit Fallback ohne Statistics Toolbox.
% X: n×d Datenmatrix, k: Nachbarzahl

X = double(X);
[n,d] = size(X);
if n <= k
    error('Für kNN-Entropie muss n>k gelten.');
end

% --- kNN-Abstände zum k-ten Nachbarn ermitteln ---
haveKnn = exist('knnsearch','file') == 2;
if haveKnn
    % Schnellweg mit Statistics Toolbox
    [~, dist] = knnsearch(X, X, 'K', k+1, 'Distance', 'euclidean');
    epsilon = dist(:, end);   % Abstand zum k-ten Nachbarn (eigener Punkt ist in Spalte 1)
else
    % Fallback: paarweise Distanzen ohne Toolbox
    % Distanzmatrix via quadratischer Form (numerisch stabilisiert)
    sq = sum(X.^2, 2);                 % n×1
    D2 = max(sq + sq' - 2*(X*X.'), 0); % n×n, quad. Distanzen (>=0)
    D  = sqrt(D2);
    D(1:n+1:end) = inf;                % Diagonale ausschließen
    D_sorted = sort(D, 2, 'ascend');
    epsilon = D_sorted(:, k);          % k-ter Nachbar (ohne Selbsttreffer)
end

% Numerische Robustheit: kein log(0)
epsilon = max(epsilon, realmin);

% Volumen der d-dimensionalen Einheitskugel: V_d = pi^(d/2) / gamma(d/2 + 1)
cd = pi^(d/2) / gamma(d/2 + 1);

% Kozachenko–Leonenko-Schätzer
H = psi(n) - psi(k) + log(cd) + d * mean(log(epsilon));
end


function D = kl_knn(P, Q, k)
% kNN-KL-Schätzer D_KL(P||Q) für d-dimensionale stetige Daten
% P: n×d, Q: m×d
P = double(P); Q = double(Q);
[n,d] = size(P);
m = size(Q,1);
if n <= k || m <= k
    error('Für kNN-KL gilt n>k und m>k.');
end

useKnnSearch = exist('knnsearch','file') == 2;
if useKnnSearch
    % Abstand zum k-ten Nachbarn in P (ohne Selbsttreffer)
    [~, distPP] = knnsearch(P, P, 'K', k+1, 'Distance', 'euclidean');
    rho = distPP(:, end);                 % n×1
    
    % Abstand zum k-ten Nachbarn aus Q
    [~, distPQ] = knnsearch(Q, P, 'K', k, 'Distance', 'euclidean');
    nu  = distPQ(:, end);                 % n×1
else
    DP = squareform(pdist(P,'euclidean')); DP(1:n+1:end) = inf;
    DQ = pdist2(P, Q, 'euclidean');
    rho = sort(DP, 2, 'ascend'); rho = rho(:, k+1);  % k+1 wg. Selbsttreffer
    nu  = sort(DQ, 2, 'ascend'); nu  = nu(:, k);     % k-ter Nachbar in Q
end

% Numerisch robust
rho = max(rho, realmin);
nu  = max(nu , realmin);

% Klassischer kNN-KL-Schätzer
D = d * mean(log(nu ./ rho)) + log(m / (n - 1));
end
