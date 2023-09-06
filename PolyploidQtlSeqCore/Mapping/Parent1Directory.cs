namespace PolyploidQtlSeqCore.Mapping
{
    /// <summary>
    /// Parent1ディレクトリ
    /// </summary>
    internal class Parent1Directory : ISampleDirectory
    {
        /// <summary>
        /// Parent1ディレクトリを作成する。
        /// </summary>
        /// <param name="dirPath">Parent1ディレクトリのPath</param>
        public Parent1Directory(string dirPath)
        {
            if (string.IsNullOrWhiteSpace(dirPath)) throw new ArgumentException(null, nameof(dirPath));

            Path = System.IO.Path.GetFullPath(dirPath);
            SampleName = System.IO.Path.GetFileName(Path);
        }

        public string Path { get; }

        public string SampleName { get; }
    }
}
