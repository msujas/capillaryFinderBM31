module to find capillaries for zscan data on BM31. Data must have headers of 'zpos i1 monitor'. Gaussian fitting is done to refine the positions.

Example usage:

```python
import capFinderBM31 as cf
#to return capillary positions and plot the results 
capPositions = cf.plotResults(<filename>)
print(capPositions)

#Or to just get the positions without plotting
capPositions = cf.getPositions(<filename>)
```

![image](https://github.com/user-attachments/assets/0f8a5ede-d6c4-45bb-ae6c-5cd30b243f60)
