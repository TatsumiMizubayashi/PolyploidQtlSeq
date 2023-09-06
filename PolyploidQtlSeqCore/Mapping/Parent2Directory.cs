namespace PolyploidQtlSeqCore.Mapping
{
    /// <summary>
    /// Parent2ディレクトリ
    /// </summary>
    internal class Parent2Directory : ISampleDirectory
    {
        /// <summary>
        /// Parent2ディレクトリを作成する。
        /// </summary>
        /// <param name="dirPath">Parent2ディレクトリのPath</param>
        public Parent2Directory(string dirPath)
        {
            if (string.IsNullOrWhiteSpace(dirPath)) throw new ArgumentException(null, nameof(dirPath));

            Path = System.IO.Path.GetFullPath(dirPath);
            SampleName = System.IO.Path.GetFileName(Path);
        }


        public string Path { get; }

        public string SampleName { get; }

    }
}
