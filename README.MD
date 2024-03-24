# Jupyter内核安装及知识库创建

- Jupyter内核
    - 参考
        - [ServiceStack Documentation](https://docs.servicestack.net/)
        - [C# Jupyter Notebooks (servicestack.net)](https://docs.servicestack.net/jupyter-notebooks-csharp#generate-c-jupyter-notebooks)
    - 查看内核列表： `jupyter kernelspec list`
    - 安装C#内核：
        - `dotnet tool install -g Microsoft.dotnet-interactive`
        - `dotnet interactive jupyter install`
    - Matlab
    账号：20229815878@cqu.edu.cn 密码：RSgpeugXWGs8egX7eR
    [https://docs.qq.com/doc/DWGJjSEdkY1lOdVJh《官网下载教程》](https://docs.qq.com/doc/DWGJjSEdkY1lOdVJh%E3%80%8A%E5%AE%98%E7%BD%91%E4%B8%8B%E8%BD%BD%E6%95%99%E7%A8%8B%E3%80%8B)
    - 安装Matlab内核：
        - `pip install jupyter-matlab-proxy`
    - 安装TypeScript内核：
        - `npm install -g tslab`
        - `tslab install --version`
        - `tslab install [--python=python3]`
    - 安装C++内核
        - 参考：[https://blog.jupyter.org/interactive-workflows-for-c-with-jupyter-fe9b54227d92](https://blog.jupyter.org/interactive-workflows-for-c-with-jupyter-fe9b54227d92)
    - 安装Go内核
    - 安装Rust内核
    - 安装Julia内核
- Jupyter创建知识库
    - 安装Jupyter Book：`pip install jupyter-book`
    - 安装ghp-import： `pip install ghp-import`
    - `jupyter-book create mynewbook/`
    - 关闭notebook执行选项，在_config.yml文件中，修改配置项：
        ```YAML
        execute:
          execute_notebooks: off
        ```
- 使用
    - `open-note`： 打开jupyter notebook
    - `open-lab`： 打开jupyter lab
    - `build <directory name>`: 构建book,生成文档复制到dist目录
    - `publish`: 发布到Github仓库的gh-pages分支，可以通过<user name>.github.io/<repo name>/<book name>访问
