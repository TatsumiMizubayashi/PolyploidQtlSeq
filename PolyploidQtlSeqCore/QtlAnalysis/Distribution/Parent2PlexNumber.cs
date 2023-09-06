namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// 親2のPlex数
    /// </summary>
    internal class Parent2PlexNumber
    {
        /// <summary>
        /// plex数の最小値
        /// </summary>
        private const int MINIMUM = 1;


        /// <summary>
        /// plex数の仮の最大値
        /// </summary>
        private const int MAXIMUM = 20;


        /// <summary>
        /// 親2のPlex数を作成する。
        /// </summary>
        /// <param name="plex">plex数</param>
        public Parent2PlexNumber(int plex)
        {
            if (plex < MINIMUM || plex > MAXIMUM) throw new ArgumentOutOfRangeException(nameof(plex));

            Value = plex;
        }

        /// <summary>
        /// Plex数を取得する。
        /// </summary>
        internal int Value { get; }

    }
}
