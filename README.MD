# Jupyter内核安装及知识库创建

知识库

下面的库和文档演示了：

- 通过Jupyter Book创建知识库
- 通过Jupyter Notebook或Jupyter Lab创建和修改记事本
- 让Jupyter记事本支持C#、Matlab、TypeScript、JavaScript等语言内核

使用Jupyter记事本创建知识库的收益包括：

- 通过记事本演示算法、模型
- 通过数据列表和数据可视化的方式跟客户讨论及验证算法和模型
- 沉淀、记录和升级知识储备
- 帮助团队成员之间、成员和客户之间的沟通、学习和交流

参考：

- Github Repo: https://github.com/richardcuick/kb-demo
- 安装文档：https://github.com/richardcuick/kb-demo/blob/main/README.MD
- 演示：https://richardcuick.github.io/kb-demo/book1/csharp-sample.html

建议：

- 大家为手中开发的程序和软件创建知识库
- 积极使用，遇到问题及时沟通解决，帮助团队、实验室迅速建立知识库
- 知识库是帮团队沉淀业务知识最好的手段，后期可以起到减少知识负载加速开发的目的。

---
安装步骤

- 安装Anaconda: [https://www.anaconda.com/download](https://www.anaconda.com/
download)
    - 使用Anaconda Prompt(Anaconda提示行，已安装Python环境)
- 创建Jupter Book
    - `jupyter-book create <book name>/`
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
        - `install -c conda-forge xeus-cling`
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