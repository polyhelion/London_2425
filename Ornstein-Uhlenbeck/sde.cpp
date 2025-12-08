#include <iostream>
#include <random>
#include <vector>
#include <optional>
#include <cmath>
#include <utility>
#include <algorithm>


#include <TApplication.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TRandom3.h>


template<typename T>
std::vector<T> arange(const T start, const T stop, const T step = T{1}) {

    std::vector<T> values;
    
    if (step == 0) {
        throw std::invalid_argument("Step size cannot be zero.");
    }

    if (step > 0) {

    	const T limit = stop - step;
    	if (limit <= 0) {
        	throw std::invalid_argument("Step size cannot be equal or bigger than the interval.");
    	}

    	for (T value = start; value < limit; value += step) {
        	values.push_back(value);
    	}

    } else {

    	const T limit = stop + step;
    	if (limit >= start) {
        	throw std::invalid_argument("Step size cannot be equal or bigger than the interval.");
    	}

    	for (T value = start; value > limit; value += step) {
        	values.push_back(value);
    	}
    }
    
    return values;
}

template <typename T, typename U>
std::vector<std::pair<T, U>> zip(const std::vector<T>& x, const std::vector<U>& y) {

    std::vector<std::pair<T, U>> result;
    std::size_t min_size = std::min(x.size(), y.size());
    result.reserve(min_size);

    std::transform(x.begin(), x.begin() + min_size, y.begin(), std::back_inserter(result),
                   [](const T& a, const U& b) { return std::make_pair(a, b); });

    return result;
}


using Time = double;


template < typename T /* floating point type data */ >
class OUProcess {

	//  Model of an Ornstein-Uhlenbeck process 

	//  dX_t = theta (mu - X_t) dt + sigma dW_t 
	//  X_0  = x0
	
	//  theta > 0  :  rate of mean reversion
	//  mu         :  long term mean
    //  sigma > 0  :  volatility
    //  (W_t)      :  a standard Brownian motion/Wiener process
    //  x0         :  initial value


	//  X_t = mu + (x0 - mu) exp(-theta t) + sigma \int_0^t exp(-theta(t-s)) dW_s

	//  E[X_t|X_0 = x0] = x0 exp(-theta t) + mu (1 - exp(-theta t))

	//  Cov(X_s, X_t) = (sigma^2 / 2 theta) [ exp(-theta |t-s|) - exp(-theta (t+s)) ]


public:

	using rd_rt = std::random_device::result_type;

	OUProcess(T theta, T mu, T sigma, T x0,
		rd_rt seed = (std::random_device())()) :
		_theta(theta), _mu(mu), _sigma(sigma),_x0(x0),
		_seed(seed),
		gen(std::mt19937(_seed)) { /* some checks go here */ }

	inline auto reversion() -> T const { return _theta; } 
	inline auto mean() -> T const { return _mu; }
	inline auto volatility() -> T const { return _sigma; }
	inline auto x0() -> T const { return _x0; }
	inline auto seed() -> T const { return _seed; }


	inline auto expectation(T x0, Time dt) const {
        return _mu + (x0 - _mu) * std::exp(-_theta * dt);
    }

    // The OU process is a continuous-time analogue of the discrete-
    // time AR(1) process. 
	// Sampling of the continuous-time solution X_t produces a Markov chain, a
	// discrete-time AR(1) process (X_k) to be precise.

	// X_{t + dt} = mu + (X_t - mu) exp(-theta dt) + 
	//                 sqrt[ (sigma^2 / (2 theta))(1 - exp(-2 theta dt))] Z
	// Z ~ N(0, 1)
    //
    // X_k = X_{k dt}


	auto sample(const Time t_0 = T{0.0}, const Time t_n = T{1.0}, 
		        const std::size_t size = 1000) 
		-> std::pair< std::vector<Time>, std::vector<T> > const {


		if (t_0 >= t_n) {
			throw std::invalid_argument("Invalid interval.");
		}
		if (size < 1) {
			throw std::invalid_argument("Invalid sample size.");
		}

	  	constexpr T zero = static_cast<T>(0.0);
    	constexpr T one  = static_cast<T>(1.0);
    	constexpr T two  = static_cast<T>(2.0);

    	const Time delta = (t_n - t_0);
    	const T scdelta  = static_cast<T>(delta);
		const Time dt = static_cast<Time>((t_n - t_0) / size);
		const T scdt = static_cast<T>(dt);

    	auto ts = arange<Time>(t_0, t_n + dt, dt);

    	const T variance = ((_sigma * _sigma) / (one * _theta)) * \
    		           (one - std::exp(-_theta * two * scdt));

    	std::normal_distribution<T> dist(zero, one);

    	// Time t;

    	std::vector<T> xs(ts.size(), zero);
        if(t_0 == zero) {
        	xs[0] = _x0;
        } else {
        	const T variance0 = ((_sigma * _sigma) / (one * _theta)) * \
    		           (one - std::exp(-_theta * two * scdelta));
        	xs[0] = _mu + (_x0 - _mu) * std::exp(-_theta * scdelta) + \
			           std::sqrt(variance0) * dist(gen);
        }

        for(std::size_t i{1}; i < ts.size(); i++) {
			// t = t_0 + (i - 1) * dt;
			xs[i] = _mu + (xs[i-1] - _mu) * std::exp(-_theta * scdt) + \
			           std::sqrt(variance) * dist(gen);
        }

    	return std::make_pair(ts, xs);

    }


	// Since Ornstein-Uhlenbeck processes have an exact discretization, other
	// discretization methods are rarely used. But for the sake of completion :

    // Euler Maruyama and Milstein method each produce a discrete-time AR(1) 
    // process which approximates the AR(1) process (X_k). 

	// X_{t + dt} = mu + (X_t - mu) dt + sigma sqrt(dt) Z
	// Z ~ N(0, 1)

	// If t_n - t_0 is too big, this implementation is problematic.
	// We need to introduce substeps. 


	auto EulerMaruyama(const Time t_0 = T{0.0}, const Time t_n = T{1.0}, 
		               const std::size_t size = 1000) 
		-> std::pair< std::vector<Time>, std::vector<T> > const {

		if (t_0 >= t_n) {
			throw std::invalid_argument("Invalid interval.");
		}
		if (size < 1) {
			throw std::invalid_argument("Invalid sample size.");
		}

		constexpr T zero = static_cast<T>(0.0);
    	constexpr T one  = static_cast<T>(1.0);

		const Time delta = (t_n - t_0);
    	const T scdelta  = static_cast<T>(delta);
    	const Time dt = static_cast<Time>((t_n - t_0) / size);
    	const T scdt = static_cast<T>(dt);

    	auto ts = arange<Time>(t_0, t_n + dt, dt);
  	
    	std::normal_distribution<T> dist(zero, one);
    	
        // Time t;

        std::vector<T> xs(ts.size(), zero);
        if(t_0 == zero) {
        	xs[0] = _x0;
        } else {
        	xs[0] = _x0 + (_theta * (_mu - _x0)) * scdelta + \
        	          _sigma * std::sqrt(scdelta) * dist(gen);
        }

        for(std::size_t i{1}; i < ts.size(); i++) {
			// t = t_0 + (i - 1) * dt;
        	xs[i] = xs[i - 1] + (_theta * (_mu - xs[i - 1])) * scdt + \
        	          _sigma * std::sqrt(scdt) * dist(gen);

        }

    	return std::make_pair(ts, xs);

	}

	auto Milstein(Time t_0, Time t_n, std::size_t size) 
		-> std::pair< std::vector<T>, std::vector<T> > const {
		// TODO
	}


	// The Runge Kutta method produces a Markov chain (but usually not an
	// AR(1) process) to approximate X_k. 


	auto RungeKutta(Time t_0, Time t_n, std::size_t size) 
		-> std::pair< std::vector<T>, std::vector<T> > const {
		// TODO
	}

	

private:

	const T _theta;
	const T _mu;
	const T _sigma;
	const T _x0;
 
	const rd_rt _seed;

	std::mt19937 gen;

};


int main(int argc, char **argv) {

	//  Run as a CERN Root application so we have the whole
	//  toolset available.
	TApplication app("Root app", &argc, argv); 

	std::random_device rd {};
	auto seed = rd();

	OUProcess<double> p(-0.5, 0.0, 0.06, 0.0);

	auto r = p.sample(0.0, 7.0, 1000);
	auto q = p.EulerMaruyama(3.0, 7.0, 1000);


    // auto results = zip(r.first, r.second);

    // for (const auto& [t, value] : results) {
    // 	std::cout << "X_t : " << t << " " << value << std::endl;
    // }

	// std::cout << p.seed() << " " << seed << '\n';


	auto c = new TCanvas("c1", "Ornstein-Uhlenbeck", 800, 600);
	auto mg = new TMultiGraph();

	//auto frame = c->DrawFrame(0.0, -4.0, 7.0, 4.0);

    auto gr1 = new TGraph(1000, r.first.data(), r.second.data());
    gr1->SetName("gr1");
    gr1->SetTitle("Sampling from the closed form transition law");
    gr1->SetLineColor(kRed);
    gr1->SetMarkerColor(kRed);
    gr1->SetLineWidth(2);
    //gr1->Draw("LP SAME");

    auto gr2 = new TGraph(1000, q.first.data(), q.second.data());
    gr2->SetName("gr2");
    gr2->SetTitle("Euler-Maruyama");
    gr2->SetLineColor(kBlue);
    gr2->SetMarkerColor(kBlue);
    gr2->SetLineWidth(2);
    //gr2->Draw("LP SAME");
    
    mg->Add(gr1, "LP");
    mg->Add(gr2, "LP");
    mg->Draw("A");

    c->Modified();
    c->Update();

    app.Run();

    return 0;

}