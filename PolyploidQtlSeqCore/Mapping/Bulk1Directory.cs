namespace PolyploidQtlSeqCore.Mapping
{
    /// <summary>
    /// Bulk1ディレクトリ
    /// </summary>
    internal class Bulk1Directory : ISampleDirectory
    {
        /// <summary>
        /// Bulk1ディレクトリを作成する。
        /// </summary>
        /// <param name="dirPath">Bulk1ディレクトリのPath</param>
        public Bulk1Directory(string dirPath)
        {
            if (string.IsNullOrWhiteSpace(dirPath)) throw new ArgumentException(null, nameof(dirPath));

            Path = System.IO.Path.GetFullPath(dirPath);
            SampleName = System.IO.Path.GetFileName(Path);
        }

        public string Path { get; }

        public string SampleName { get; }

    }
}
