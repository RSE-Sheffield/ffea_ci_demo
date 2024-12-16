## **Git branch organisation of FFEA**

#### **Ben Hanson:**

"I think an alpha, beta, gamma branch structure would be best i.e. whenever you
work on a new feature, first switch to alpha to get the most recent version, 
then create a new branch and make changes in there. When you're finished, do
some quick tests then merge back into alpha. Devs all work on alpha, branching
off and merging back in. Then bugs will appear on the alpha branch, we can just
branch off and bug fix just as with any feature. When we're ready for a release,
only then do we merge alpha into beta, and give beta to people to test. 
Potentially have a gamma test phase too for last minute patches. Then merge into
master when we're absolutely sure everything is sound and tag it. From there, 
master would just big a big list of tagged stable merges, you see?"

* **alpha**: cutting edge development
* **beta**: testing, bug fixes
* **gamma**: hotfixes
* **master**: stable releases

[Guide](https://nickskelton.medium.com/iframe-width-420-height-315-src-https-www-youtube-com-embed-fm3d7byzlic-8a81158c1110)
