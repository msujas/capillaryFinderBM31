module to find capillaries for zscan data on BM31. Data must have 3 columns and headers of 'zpos i1 monitor'. It searches for the last scan in the file, i.e. the last instance of the column header. Gaussian fitting is done to refine the positions. Install using 'pip install -e .' in the directory.

Example usage:

```python
import capFinderBM31 as cf
#to return capillary positions and plot the results 
capPositions = cf.plotResults(<filename>)
print(capPositions)

#Or to just get the positions without plotting
capPositions = cf.getPositions(<filename>)
```

There is also a python executable which can be used in the terminal. E.g.:
```
plotZscan zscan.dat
```
Use --help to display help.

![image](https://github.com/user-attachments/assets/0f8a5ede-d6c4-45bb-ae6c-5cd30b243f60)
