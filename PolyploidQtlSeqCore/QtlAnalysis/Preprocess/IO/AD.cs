namespace PolyploidQtlSeqCore.QtlAnalysis.Preprocess.IO
{
    /// <summary>
    /// AD
    /// </summary>
    internal class AD
    {
        private static readonly char[] _delimiter = [','];
        private const string NO_DATA = ".";

        /// <summary>
        /// NoData ADインスタンスを作成する。
        /// </summary>
        public AD() 
            : this(NO_DATA)
        {
        }


        /// <summary>
        /// ADインスタンスを作成する。
        /// </summary>
        /// <param name="ad">ADの値</param>
        public AD(string ad)
        {
            ArgumentException.ThrowIfNullOrEmpty(ad);

            Value = ad;
            var isNoData = ad.Contains(NO_DATA);
            var readCounts = isNoData ? new[] {0, 0} : [.. ad.Split(_delimiter).Select(x => int.Parse(x))];
            RefCount = readCounts[0];
            AltCount = readCounts[1];
            Depth = readCounts.Sum();
        }

        /// <summary>
        /// ADの値を取得する。
        /// </summary>
        public string Value { get; }

        /// <summary>
        /// Ref型アレルリード数を取得する。
        /// </summary>
        public int RefCount { get; }

        /// <summary>
        /// Alt型アレルリード数を取得する。
        /// </summary>
        public int AltCount { get; }

        /// <summary>
        /// Depthを取得する。
        /// </summary>
        public int Depth { get; }
    }
}
