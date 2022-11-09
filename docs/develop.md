
## Development
1. Start with the latest development branch:
   ```bash
   $ git checkout dev
   $ git fetch origin
   $ git reset --hard origin/dev
   ```
2. Create a branch with your feature or change:
   ```bash
   $ git checkout -b feature_name
   ```
3. Write code, add files, and commit the changes:
   ```bash
   $ git status
   $ git add <some-files>
   $ git commit -m "Added my cool new features"
   ```
4. Rebase frequently to incorporate upstream changes:
   ```bash
   git fetch origin dev
   git rebase origin/dev
   ```
5. Test your code locally and write proper documentation.
6. Once you are done, push your changes back to GitHub:
   ```bash
   $ git push -u origin feature_name
   ```
7. Create a Pull Request in GitHub between your feature branch and development
8. Your code will be reviewed.
9. Once approved, merge your branch to the Development branch, and delete the
   feature branch.
