namespace PolyploidQtlSeqCore.Mapping
{
    /// <summary>
    /// Bulk2ディレクトリ
    /// </summary>
    internal class Bulk2Directory : ISampleDirectory
    {
        /// <summary>
        /// Bulk2ディレクトリを作成する。
        /// </summary>
        /// <param name="dirPath">Bulk2ディレクトリのPath</param>
        public Bulk2Directory(string dirPath)
        {
            if (string.IsNullOrWhiteSpace(dirPath)) throw new ArgumentException(null, nameof(dirPath));

            Path = System.IO.Path.GetFullPath(dirPath);
            SampleName = System.IO.Path.GetFileName(Path);
        }


        public string Path { get; }

        public string SampleName { get; }
    }
}
