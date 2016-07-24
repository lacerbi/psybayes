# psybayes
Adaptive Bayesian stimulus placement of psychometric function for MATLAB.

`psybayes.m` implements Kontsevich and Tyler's (1999) Bayesian adaptive Ψ method for estimation of parameters of the psychometric function via maximization of information gain (including lapse, see Prins 2012). It also supports the marginal-Ψ method by Prins (2013). 

See `psytest.m` for documentation and a working usage example.

The output of **Psybayes** can be plotted with `psybayes_plot.m`:

![alt text](http://luigiacerbi.com/wp-content/uploads/2016/07/psybayes_fig.png "psybayes demo")


### Notes

This is a fairly bare implementation of Bayesian adaptive stimulus placement. For an array of adaptive procedures you could look at the [Palamedes](http://www.palamedestoolbox.org/) toolbox by Nicolaas Prins and Frederick Kingdom, which has some pretty neat [demo](http://www.palamedestoolbox.org/pal_ampm_demo.html).

### References

Kontsevich, L. L., & Tyler, C. W. (1999). "Bayesian adaptive estimation of psychometric slope and threshold". *Vision research*, 39(16), 2729-2737. ([link](http://www.sciencedirect.com/science/article/pii/S0042698998002855))

Prins, N. (2012). The adaptive psi method and the lapse rate. *Journal of Vision*, 12(9), 322-322. ([link](http://jov.arvojournals.org/article.aspx?articleid=2140969))

Prins, N. (2013). "The psi-marginal adaptive method: How to give nuisance parameters the attention they deserve (no more, no less)". *Journal of vision*, 13(7), 3-3. ([link](http://jov.arvojournals.org/article.aspx?articleid=2121598))

### License

This is free, open source software released under the GPL-3 License. See LICENSE.txt for details.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Use this software at your own risk.
